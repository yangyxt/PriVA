#!/usr/bin/env python3

import pysam
import numpy as np
from collections import defaultdict
from typing import Dict, List
from multiprocessing import Pool
import sys
import pickle

# Define the nested defaultdict factory at the module level
def nested_defaultdict():
    return defaultdict(nested_defaultdict)


def to_plain_dict(obj):
    """Strip every defaultdict out of a nested structure, keeping the data.

    WHY THE OUTPUT MUST NOT CONTAIN A defaultdict
    =============================================

    pickle does not store a function; it stores a module name and an attribute
    name, and looks the pair up on load. A defaultdict therefore pickles a
    reference to its factory -- and because this file is normally RUN rather
    than imported, that reference came out as ``__main__.nested_defaultdict``.

    On load, pickle evaluates ``getattr(sys.modules["__main__"], "nested_
    defaultdict")``, which asks whichever file is the ENTRY POINT of the reading
    program, not the file calling pickle.load. That made the caches readable
    from one entry point and not another:

        python am_pick_intolerant_domains.py --pickle_file ...      worked
          that file imports the factory at its top, and it is __main__

        python prepare_pm1_regions.py ...                           FAILED
          it imports analyze_domain_data and calls it, so __main__ is
          prepare_pm1_regions, which never names the factory

    The second is a real path in install_utils.sh, and it raised
    ``AttributeError: Can't get attribute 'nested_defaultdict'`` the moment the
    intolerant-domain pickle needed regenerating.

    Converting to plain dicts before dumping removes the reference entirely, so
    the caches are readable by any program. Nothing is lost: readers use
    .items(), .values(), .get() and ``in``, so none depends on a missing key
    being created on access.

    The four deployed protein-domain caches were converted in place on
    2026-07-30 with contents verified unchanged; this keeps rebuilds clean.
    """
    if isinstance(obj, dict):          # covers defaultdict, which subclasses dict
        return {k: to_plain_dict(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [to_plain_dict(v) for v in obj]
    if isinstance(obj, tuple):
        return tuple(to_plain_dict(v) for v in obj)
    return obj

# VEP DOMAINS sources that are structural mappings, not sequence-homology
# domain annotations. Filtering them here avoids inflating scores_dict with
# entries that span most of the protein (PDB/AFDB coordinate mappings) or
# represent low-complexity / disordered regions (MobiDB_lite).
_NON_DOMAIN_SOURCES = frozenset({
    "PDB-ENSP_mappings",
    "AFDB-ENSP_mappings",
    "MobiDB_lite",
})


class DomainAMScoreCollector:
    '''
    {
        <ENSG_ID>: {  # Key: Ensembl Gene ID (string, e.g., 'ENSG00000123456')
            <DB_NAME>: {  # Key: Domain Database Name (string, e.g., 'PANTHER', 'Pfam')
            <LEVEL_1_ID>: {  # Key: First level of domain ID (string, e.g., 'PTHR26451', 'PF00001')
                # --- This level contains transcript info and the score distribution ---
                <TRANSCRIPT_ID_1>: set(<EXON_ID_A>, <EXON_ID_B>, ...), # Key: Ensembl Transcript ID (string), Value: Set of Exon IDs (strings like '1/10')
                <TRANSCRIPT_ID_2>: set(<EXON_ID_C>, ...),
                # ... potentially more transcript IDs ...
                'distribution': np.ndarray([<score1>, <score2>, ...]), # Key: literal string 'distribution', Value: NumPy array of AlphaMissense scores (floats) for this domain level

                # --- If the domain has further hierarchy (e.g., PANTHER subfamilies) ---
                <LEVEL_2_ID>: { # Key: Second level of domain ID (string, e.g., 'SF72')
                # --- This level ALSO contains transcript info and its own score distribution ---
                <TRANSCRIPT_ID_1>: set(<EXON_ID_A>, <EXON_ID_D>, ...), # Transcript info specific to this sub-level
                <TRANSCRIPT_ID_3>: set(<EXON_ID_E>, ...),
                # ... potentially more transcript IDs ...
                'distribution': np.ndarray([<score3>, <score4>, ...]), # NumPy array of scores specific to this sub-level (LEVEL_1 + LEVEL_2)

                # --- And potentially more levels ---
                <LEVEL_3_ID>: {
                    # ... structure repeats ...
                    'distribution': np.ndarray([...])
                                }
                            }
                        }
                        # ... potentially more LEVEL_1_IDs ...
                        }
                        # ... potentially more DB_NAMEs ...
                    }
                    # ... potentially more ENSG_IDs ...
    }
    '''
    def __init__(self, vcf_path: str):
        """Initialize collector with VCF file path."""
        self.vcf_path = vcf_path
        self.scores_dict = nested_defaultdict()
        
    def _parse_domain_hierarchy(self, domains_str: str) -> List[tuple]:
        """Parse domain string into hierarchical components.
        
        Args:
            domains_str: String like 'PANTHER:PTHR26451&PANTHER:PTHR26451:SF72'
            
        Returns:
            List of tuples containing (db_name, hierarchy_list)
            e.g., [('PANTHER', ['PTHR26451']), ('PANTHER', ['PTHR26451', 'SF72'])]
        """
        if not domains_str or domains_str == '':
            return []
            
        domain_entries = domains_str.split('&')
        parsed_domains = []
        
        for entry in domain_entries:
            parts = entry.split(':')
            db_name = parts[0]
            if db_name in _NON_DOMAIN_SOURCES:
                continue

            # Build hierarchical structure
            for i in range(1, len(parts)):
                hierarchy = parts[1:i+1]  # Include all levels up to current
                parsed_domains.append((db_name, hierarchy))
                
        return parsed_domains

    def collect_scores(self, region=None):
        """Process VCF file and collect AlphaMissense pathogenicity scores.
        
        Args:
            region: Optional chromosome/region string for region-based fetching.
        """
        vcf = pysam.VariantFile(self.vcf_path)
        
        # Get CSQ format from header
        csq_format = str(vcf.header).split('Format: ')[-1].strip('>"').split('|')
        domains_idx = csq_format.index('DOMAINS')
        gene_idx = csq_format.index('Gene')
        exon_idx = csq_format.index('EXON')
        protein_pos_idx = csq_format.index('Protein_position')
        feature_type_idx = csq_format.index('Feature_type')
        transcript_idx = csq_format.index('Feature')
        
        records = vcf.fetch(region) if region else vcf
        for record in records:
            try:
                am_score = float(record.info['AM_PATHOGENICITY'])
                try:
                    gnomAD_AF = float(record.info['AF_grpmax_joint'][0])
                except (KeyError, ValueError, IndexError):
                    gnomAD_AF = 0.0
                prot_var = record.info['PVAR']
                # Extract amino acid position from PVAR (assuming format like 'p.Arg123Ser')
                aa_pos = prot_var[:-1]

                csq_value = record.info['CSQ']
                if isinstance(csq_value, tuple):
                    transcript_annotations = csq_value
                else:
                    transcript_annotations = csq_value.split(',')

                for csq_entry in transcript_annotations:
                    fields = csq_entry.split('|')
                    
                    # Skip if not a transcript or not an ENST transcript
                    if (fields[feature_type_idx] != 'Transcript' or 
                        not fields[transcript_idx].startswith('ENST')):
                        continue
                    
                    domains_str = fields[domains_idx]
                    ensg_id = fields[gene_idx]
                    transcript_id = fields[transcript_idx]
                    exon_id = fields[exon_idx]
                    protein_pos_str = fields[protein_pos_idx]
                    
                    # Parse integer protein position from VEP format "150/500" or "150-151/500"
                    aa_pos_int = None
                    if protein_pos_str:
                        try:
                            aa_pos_int = int(protein_pos_str.split('/')[0].split('-')[0])
                        except ValueError:
                            pass
                    
                    if not domains_str or not ensg_id:
                        continue
                    
                    # Process each domain annotation
                    for db_name, hierarchy in self._parse_domain_hierarchy(domains_str):
                        current_dict = self.scores_dict[ensg_id][db_name]
                        
                        # Navigate through domain hierarchy
                        for level in hierarchy:
                            current_dict = current_dict[level]
                            if transcript_id not in current_dict:
                                current_dict[transcript_id] = {'exons': set(), 'aa_min': float('inf'), 'aa_max': 0}
                            current_dict[transcript_id]['exons'].add(exon_id)
                            if aa_pos_int is not None:
                                if aa_pos_int < current_dict[transcript_id]['aa_min']:
                                    current_dict[transcript_id]['aa_min'] = aa_pos_int
                                if aa_pos_int > current_dict[transcript_id]['aa_max']:
                                    current_dict[transcript_id]['aa_max'] = aa_pos_int
                            # Initialize distribution if not exists
                            if 'distribution' not in current_dict:
                                current_dict['distribution'] = {}
                            # Add score to distribution
                            if prot_var not in current_dict['distribution']:
                                current_dict['distribution'][prot_var] = am_score

                            if 'gnomAD_distribution' not in current_dict:
                                current_dict['gnomAD_distribution'] = {}
                            if prot_var not in current_dict['gnomAD_distribution']:
                                current_dict['gnomAD_distribution'][prot_var] = gnomAD_AF
                            elif gnomAD_AF > current_dict['gnomAD_distribution'][prot_var]:
                                current_dict['gnomAD_distribution'][prot_var] = gnomAD_AF

                            if 'gnomAD_residue_distribution' not in current_dict:
                                current_dict['gnomAD_residue_distribution'] = {}
                            if aa_pos not in current_dict['gnomAD_residue_distribution']:
                                current_dict['gnomAD_residue_distribution'][aa_pos] = [gnomAD_AF]
                            else:
                                current_dict['gnomAD_residue_distribution'][aa_pos].append(gnomAD_AF)

                            if 'min_distribution' not in current_dict:
                                current_dict['min_distribution'] = {}
                            if aa_pos not in current_dict['min_distribution']:
                                current_dict['min_distribution'][aa_pos] = am_score
                            elif am_score < current_dict['min_distribution'][aa_pos]:
                                # Keep the lowest score for this AA position
                                current_dict['min_distribution'][aa_pos] = am_score

                            if 'max_distribution' not in current_dict:
                                current_dict['max_distribution'] = {}
                            if aa_pos not in current_dict['max_distribution']:
                                current_dict['max_distribution'][aa_pos] = am_score
                            elif am_score > current_dict['max_distribution'][aa_pos]:
                                # Keep the highest score for this AA position
                                current_dict['max_distribution'][aa_pos] = am_score
            except (KeyError, ValueError) as e:
                print(f"Error processing record {record.chrom}:{record.pos} - {str(e)}")
                continue

    def finalize_scores(self):
        """Convert score dict to numpy arrays."""
        def convert_nested(d):
            for k, v in d.items():
                if k == 'distribution' and isinstance(v, dict):
                    d[k] = np.array(list(v.values()))
                elif k == 'min_distribution' and isinstance(v, dict):
                    d[k] = np.array(list(v.values()))
                elif k == 'max_distribution' and isinstance(v, dict):
                    d[k] = np.array(list(v.values()))
                elif k == 'gnomAD_distribution' and isinstance(v, dict):
                    d[k] = np.array(list(v.values()))
                elif isinstance(v, dict):
                    convert_nested(v)
        
        convert_nested(self.scores_dict)
    
    @staticmethod
    def load_scores_from_pickle(pickle_path: str) -> Dict:
        """Load scores from the pickle file."""
        with open(pickle_path, 'rb') as f:
            data = pickle.load(f)
        return data


def _collect_worker(args):
    """Multiprocessing worker: collect scores for a single chromosome."""
    vcf_path, chrom = args
    collector = DomainAMScoreCollector(vcf_path)
    collector.collect_scores(region=chrom)
    return dict(collector.scores_dict)


def _merge_domain_level(target, source):
    """Recursively merge source domain-level dict into target in place."""
    for key, value in source.items():
        if key not in target:
            target[key] = value
            continue
        if isinstance(key, str) and key.startswith('ENST'):
            # Transcript entry — merge exon sets and aa ranges
            tv = target[key]
            if isinstance(value, dict) and 'exons' in value:
                tv['exons'].update(value['exons'])
                tv['aa_min'] = min(tv['aa_min'], value['aa_min'])
                tv['aa_max'] = max(tv['aa_max'], value['aa_max'])
            elif isinstance(value, set):
                tv.update(value)
        elif key == 'distribution':
            for k, v in value.items():
                if k not in target[key]:
                    target[key][k] = v
        elif key == 'gnomAD_distribution':
            for k, v in value.items():
                if k not in target[key] or v > target[key][k]:
                    target[key][k] = v
        elif key == 'gnomAD_residue_distribution':
            for k, v in value.items():
                if k not in target[key]:
                    target[key][k] = list(v)
                else:
                    target[key][k].extend(v)
        elif key == 'min_distribution':
            for k, v in value.items():
                if k not in target[key] or v < target[key][k]:
                    target[key][k] = v
        elif key == 'max_distribution':
            for k, v in value.items():
                if k not in target[key] or v > target[key][k]:
                    target[key][k] = v
        elif isinstance(value, dict):
            _merge_domain_level(target[key], value)


def _merge_scores(target, source):
    """Merge a per-chromosome scores_dict into the combined target."""
    for gene_id, gene_data in source.items():
        if gene_id not in target:
            target[gene_id] = gene_data
        else:
            for db_name, domain_data in gene_data.items():
                if db_name not in target[gene_id]:
                    target[gene_id][db_name] = domain_data
                else:
                    _merge_domain_level(target[gene_id][db_name], domain_data)


def _finalize_scores_dict(scores_dict):
    """Convert distribution dicts to numpy arrays (same as DomainAMScoreCollector.finalize_scores)."""
    def convert_nested(d):
        for k, v in d.items():
            if k == 'distribution' and isinstance(v, dict):
                d[k] = np.array(list(v.values()))
            elif k == 'min_distribution' and isinstance(v, dict):
                d[k] = np.array(list(v.values()))
            elif k == 'max_distribution' and isinstance(v, dict):
                d[k] = np.array(list(v.values()))
            elif k == 'gnomAD_distribution' and isinstance(v, dict):
                d[k] = np.array(list(v.values()))
            elif k == 'gnomAD_residue_distribution' and isinstance(v, dict):
                pass  # Keep as dict {aa_pos: [AF1, AF2, ...]} — do NOT flatten
            elif isinstance(v, dict):
                convert_nested(v)
    convert_nested(scores_dict)


def main(vcf_path, output_pickle, threads=1):
    """
    Notice that the input VCF file must be the AlphaMissense VCF file annotated by VEP with the --domains argument.
    """
    if threads > 1:
        vcf = pysam.VariantFile(vcf_path)
        chroms = [c for c in vcf.header.contigs]
        vcf.close()
        n_workers = min(threads, len(chroms))
        print(f"Processing {len(chroms)} chromosomes with {n_workers} workers", file=sys.stderr)
        
        with Pool(n_workers) as pool:
            results = pool.map(_collect_worker, [(vcf_path, c) for c in chroms])
        
        merged = {}
        for result in results:
            _merge_scores(merged, result)
        
        _finalize_scores_dict(merged)
        scores_dict = merged
    else:
        collector = DomainAMScoreCollector(vcf_path)
        collector.collect_scores()
        collector.finalize_scores()
        scores_dict = collector.scores_dict
    
    if output_pickle:
        with open(output_pickle, 'wb') as f:
            pickle.dump(to_plain_dict(scores_dict), f,
                        protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Saved scores to {output_pickle}", file=sys.stderr)

    return scores_dict
        

if __name__ == "__main__":
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    main(sys.argv[1], sys.argv[2], threads)