import xml.etree.ElementTree as ET
import gzip
import sys
import time
import re
import numpy as np
import pysam # Assuming pysam is available
from typing import Dict, Tuple, List, Optional, Set, Any
import pickle
import logging
from collections import defaultdict
import json
import os
import csv
self_directory = os.path.dirname(os.path.abspath(__file__))

# --- Logger Setup ---
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)
# Check if handlers already exist to prevent duplication
if not logger.handlers:
    console_handler = logging.StreamHandler(sys.stderr) # Log to stderr
    console_handler.setLevel(logging.ERROR)
    formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(name)s:%(funcName)s:%(lineno)s:%(message)s")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
# --- End Logger Setup ---

# --- Constants ---
MAX_EXAMPLES_PER_DB = 100 # Limit memory usage for examples

# --- Type Hint Definitions for Clarity ---
GOTermInfo = Tuple[str, str, str]
InterProMapValue = Tuple[str, Optional[str], Optional[str], str, List[GOTermInfo], str]
FinalMapping = Dict[Tuple[str, str], List[InterProMapValue]]
# --- End Type Hint Definitions ---

# --- Class Definition: DomainNormalizer ---
class DomainNormalizer():
    """
    Analyzes domain annotations from VEP VCF files and InterPro XML data
    to facilitate the creation of curated mappings for database names and
    accession ID formats between the two sources.

    Workflow (within a single script run as requested):
    1. Instantiate the class.
    2. Call parse_vcf_file() to collect VEP formats.
    3. Call parse_interpro_xml_formats() to collect InterPro formats.
    4. Call generate_comparison_report() to output the collected info for review.
    5. **Note:** Normalization methods require manual curation based on the report.
    """
    def __init__(self, logger = logger):
        """Initialize the normalizer."""
        self.vep_db_names: Set[str] = set()
        self.vep_db_id_examples: Dict[str, Set[str]] = defaultdict(set)
        self.vep_db_counts: Dict[str, int] = defaultdict(int)
        self.interpro_db_names: Set[str] = set()
        self.interpro_db_id_examples: Dict[str, Set[str]] = defaultdict(set)
        self.interpro_db_counts: Dict[str, int] = defaultdict(int)
        self._vep_to_ipr_map: Dict[str, str] = {}
        self._ipr_to_vep_map: Dict[str, str] = {}
        self._vep_to_ipr_id_rules: Dict = {}
        self._ipr_to_vep_id_rules: Dict = {}
        self.logger = logger
        self.load_curated_mappings()
        self.load_curated_functional_domains()
        self.logger.info("DomainNormalizer initialized.")

    def parse_vcf_file(self, vcf_path: str) -> None:
        """ Parses VCF to collect VEP domain formats. """
        self.logger.info(f"Starting VEP domain format collection from: {vcf_path}")
        record_count = 0
        processed_domains_total = 0
        start_time = time.time()
        try:
            vcf = pysam.VariantFile(vcf_path, "r")
            csq_header_info = vcf.header.info.get('CSQ')
            if csq_header_info is None or not csq_header_info.description:
                logger.error("CSQ INFO field description not found in VCF header.")
                return
            csq_header_desc = csq_header_info.description
            if 'Format: ' not in csq_header_desc:
                logger.error("Could not parse CSQ format from header description.")
                return
            field_names = csq_header_desc.split('Format: ')[1].strip('"').split('|')
            csq_fields = {name: i for i, name in enumerate(field_names)}
            if 'DOMAINS' not in csq_fields:
                logger.error("DOMAINS field not found in CSQ format description.")
                return
            domains_idx = csq_fields['DOMAINS']

            for record in vcf:
                record_count += 1
                try:
                    if 'CSQ' in record.info:
                        csq_entries = record.info['CSQ']
                        if not isinstance(csq_entries, tuple):
                            csq_entries = (csq_entries,)
                        for csq_entry in csq_entries:
                             if isinstance(csq_entry, str):
                                fields = csq_entry.split('|')
                                if len(fields) > domains_idx and fields[domains_idx]:
                                    domains_str = fields[domains_idx]
                                    processed_domains_total += self._parse_domains_vep(domains_str)
                             else:
                                self.logger.warning(f"Skipping non-string CSQ entry in record {record.id or record_count}")
                except Exception as record_err:
                    self.logger.warning(f"Error processing record #{record_count} (ID: {record.id}): {record_err}")
                if record_count % 500000 == 0:
                    elapsed = time.time() - start_time
                    self.logger.info(f"Processed {record_count} VCF records... ({elapsed:.1f}s)")

            vcf.close()
            elapsed = time.time() - start_time
            self.logger.info(f"Finished VCF processing. Processed {record_count} records in {elapsed:.2f} seconds.")
            self.logger.info(f"Found {len(self.vep_db_names)} unique VEP domain DB names.")
            self.logger.info(f"Collected {processed_domains_total} total domain annotations.")
        except ImportError:
            self.logger.error("pysam library not found. Please install pysam to parse VCF files.")
            raise
        except FileNotFoundError:
            self.logger.error(f"VCF file not found at: {vcf_path}")
            raise
        except Exception as e:
            self.logger.exception(f"An unexpected error occurred parsing VCF file: {vcf_path}")
            raise

    def _parse_domains_vep(self, domains_str: str) -> int:
        """ Parses domain string from VEP annotation and updates VEP stats. """
        count = 0
        if not domains_str: return count
        try:
            domains = domains_str.split('&')
            for domain in domains:
                parts = domain.split(':', 1)
                if len(parts) == 2:
                    db_name, domain_id = parts
                    if db_name and domain_id:
                        self.vep_db_names.add(db_name)
                        if len(self.vep_db_id_examples[db_name]) < MAX_EXAMPLES_PER_DB:
                            self.vep_db_id_examples[db_name].add(domain_id)
                        self.vep_db_counts[db_name] += 1
                        count += 1
        except Exception as e:
            self.logger.warning(f"Error parsing VEP domain substring '{domains_str[:50]}...': {e}")
        return count

    def parse_interpro_xml_formats(self, xml_file_path: str) -> None:
        """ Parses InterPro XML to collect member DB names and example keys. """
        self.logger.info(f"Starting InterPro format collection from: {xml_file_path}")
        interpro_entry_count = 0
        processed_signatures_total = 0
        start_time = time.time()
        try:
            with gzip.open(xml_file_path, 'rb') as f:
                context = ET.iterparse(f, events=('end',))
                _, root = next(context)
                for event, elem in context:
                    if elem.tag == 'interpro':
                        interpro_entry_count += 1
                        try:
                            member_list = elem.find('member_list')
                            if member_list is not None:
                                for db_xref in member_list.findall('db_xref'):
                                    db_name = db_xref.get('db')
                                    db_key = db_xref.get('dbkey')
                                    if db_name and db_key:
                                        processed_signatures_total += 1
                                        self.interpro_db_names.add(db_name)
                                        if len(self.interpro_db_id_examples[db_name]) < MAX_EXAMPLES_PER_DB:
                                            self.interpro_db_id_examples[db_name].add(db_key)
                                        self.interpro_db_counts[db_name] += 1
                        except Exception as inner_e:
                            ipr_id = elem.get('id', 'N/A')
                            self.logger.warning(f"Error processing members for InterPro entry {ipr_id}: {inner_e}")
                        finally:
                             elem.clear()
                        if interpro_entry_count % 50000 == 0:
                            elapsed = time.time() - start_time
                            self.logger.info(f"Processed {interpro_entry_count} InterPro entries (format scan)... ({elapsed:.1f}s)")
            elapsed = time.time() - start_time
            self.logger.info(f"Finished InterPro XML format collection. Processed {interpro_entry_count} entries in {elapsed:.2f} seconds.")
            self.logger.info(f"Found {len(self.interpro_db_names)} unique InterPro member DB names.")
            self.logger.info(f"Collected {processed_signatures_total} total member signature annotations.")
        except FileNotFoundError:
            self.logger.error(f"Error: File not found at {xml_file_path}")
            raise
        except ET.ParseError as e:
            self.logger.exception(f"Error: Failed to parse XML file")
            raise
        except Exception as e:
            self.logger.exception(f"An unexpected error occurred parsing InterPro XML for formats")
            raise

    def generate_comparison_report(self, output_path: str) -> None:
        """ Generates comparison report for manual curation. """
        self.logger.info(f"Generating comparison report to: {output_path}")
        try:
            with open(output_path, 'w') as f:
                # (Content of report generation as in previous version)
                f.write("VEP vs InterPro Domain Database Format Comparison Report\n")
                f.write("=======================================================\n\n")
                f.write("Purpose: Review the database names and example IDs from both VEP annotations\n")
                f.write("         and the InterPro XML file (<member_list>) to establish curated\n")
                f.write("         mappings for both database names and necessary ID format adjustments.\n\n")
                f.write("--- VEP Annotations Found ---\n")
                f.write(f"Total unique VEP DB names: {len(self.vep_db_names)}\n\n")
                if not self.vep_db_names: f.write("(No VEP domain data collected)\n\n")
                else:
                     sorted_vep_dbs = sorted(list(self.vep_db_names))
                     for db in sorted_vep_dbs:
                         f.write(f"VEP DB Name: {db}\n")
                         f.write(f"  Occurrences: {self.vep_db_counts.get(db, 0)}\n")
                         examples = list(self.vep_db_id_examples[db])[:5]
                         f.write(f"  Example IDs: {', '.join(examples)}\n\n")
                f.write("\n--- InterPro Member Databases Found (<member_list>) ---\n")
                f.write(f"Total unique InterPro member DB names: {len(self.interpro_db_names)}\n\n")
                if not self.interpro_db_names: f.write("(No InterPro member data collected)\n\n")
                else:
                     sorted_ipr_dbs = sorted(list(self.interpro_db_names))
                     for db in sorted_ipr_dbs:
                         f.write(f"InterPro DB Name ('db' attribute): {db}\n")
                         f.write(f"  Occurrences: {self.interpro_db_counts.get(db, 0)}\n")
                         examples = list(self.interpro_db_id_examples[db])[:5]
                         f.write(f"  Example IDs ('dbkey' attribute): {', '.join(examples)}\n\n")
                f.write("\n--- Curation Action Required ---\n")
                f.write("1. Review VEP names vs InterPro names and create VEP -> InterPro DB mapping.\n")
                f.write("2. Review VEP IDs vs InterPro IDs and define ID normalization rules (e.g., strip :SFxxx from PANTHER).\n")
                f.write("3. Implement mappings and rules in DomainNormalizer methods.\n")

            self.logger.info(f"Comparison report saved successfully to {output_path}")
        except Exception as e:
            self.logger.exception(f"Failed to generate or save comparison report to {output_path}")

    def load_curated_mappings(self, json_file_path = os.path.join(os.path.dirname(self_directory), 'data', 'InterPro', 'InterPro_domain_norm_map.json')) -> None:
        """ Loads the manually curated mappings from a JSON file. """
        try:
            with open(json_file_path, 'r') as f:
                mapping_data = json.load(f)
                
            self._vep_to_ipr_map = mapping_data['database_name_mappings']['vep_to_interpro']
            self._ipr_to_vep_map = mapping_data['database_name_mappings']['interpro_to_vep']
            self._vep_to_ipr_id_rules = mapping_data['accession_id_normalization']['vep_to_interpro_key']
            self._ipr_to_vep_id_rules = mapping_data['accession_id_normalization']['interpro_key_to_vep_id']
            
            self.logger.info(f"Loaded domain mappings from {json_file_path}")
            self.logger.info(f"Loaded {len(self._vep_to_ipr_map)} VEP to InterPro DB mappings")
            self.logger.info(f"Loaded {len(self._ipr_to_vep_map)} InterPro to VEP DB mappings")
        except Exception as e:
            self.logger.error(f"Error loading domain mappings from {json_file_path}: {e}")
            raise

    def normalize_vep_to_interpro(self, vep_db_name: str, vep_id: str) -> Optional[Tuple[str, str]]:
        """ Normalizes VEP DB name and ID to InterPro format. """
        if not self._vep_to_ipr_map or not self._vep_to_ipr_id_rules:
            self.logger.error("Curated mappings not loaded. Call load_curated_mappings() first.")
            return None
            
        # Map database name
        interpro_db = self._vep_to_ipr_map.get(vep_db_name)
        if not interpro_db:
            return None  # No mapping for this database
            
        # Apply ID normalization rules
        rule = self._vep_to_ipr_id_rules.get(vep_db_name, self._vep_to_ipr_id_rules['__default__'])
        rule_type = rule.get('rule_type')
        
        normalized_id = vep_id  # Default to original ID
        
        if rule_type == 'identity':
            # Use ID as is
            pass
            
        elif rule_type == 'regex_capture':
            # Extract portion of ID using regex
            pattern = rule.get('pattern')
            group = rule.get('group', 1)
            if pattern:
                match = re.match(pattern, vep_id)
                if match and len(match.groups()) >= group:
                    normalized_id = match.group(group)
                else:
                    self.logger.debug(f"Regex pattern {pattern} did not match VEP ID {vep_id}")
                    return None
            
        elif rule_type == 'prefix_add':
            # Add prefix to ID
            prefix = rule.get('prefix', '')
            normalized_id = f"{prefix}{vep_id}"
            
        elif rule_type == 'prefix_check':
            # Check if prefix exists, add if missing
            prefix = rule.get('prefix', '')
            action = rule.get('action_if_missing')
            if not vep_id.startswith(prefix) and action == 'prepend':
                normalized_id = f"{prefix}{vep_id}"
        
        if not normalized_id:
            return None
        
        self.logger.debug(f"Normalized VEP anno domain {vep_db_name}:{vep_id} to InterPro XML domain {interpro_db}:{normalized_id}")
        return (interpro_db, normalized_id)


    def normalize_interpro_to_vep(self, interpro_db_name: str, interpro_dbkey: str) -> Optional[Tuple[str, str]]:
        """ Normalizes InterPro DB name and key back to VEP format. """
        if not self._ipr_to_vep_map or not self._ipr_to_vep_id_rules:
            self.logger.error("Curated mappings not loaded. Call load_curated_mappings() first.")
            return None
            
        # Map database name
        vep_db_name = self._ipr_to_vep_map.get(interpro_db_name.upper())
        if not vep_db_name:
            return None  # No mapping for this database
            
        # Apply ID normalization rules
        rule = self._ipr_to_vep_id_rules.get(interpro_db_name.upper(), self._ipr_to_vep_id_rules['__default__'])
        rule_type = rule.get('rule_type')
        
        normalized_id = interpro_dbkey  # Default to original ID
        
        if rule_type == 'identity':
            # Use ID as is
            pass
            
        elif rule_type == 'regex_substitute':
            # Replace part of ID using regex
            pattern = rule.get('pattern')
            replacement = rule.get('replacement', '')
            if pattern:
                normalized_id = re.sub(pattern, replacement, interpro_dbkey)
        
        if not normalized_id:
            return None
            
        return (vep_db_name, normalized_id)

    

    def query_interpro_entry_vep_anno(self, vep_anno_str: str, interpro_entry_map_dict: FinalMapping, vep_pred_domains = { "Cleavage_site_(Signalp)",
                                                                                                              "Coiled-coils_(Ncoils)",
                                                                                                              "Low_complexity_(Seg)",
                                                                                                              "Transmembrane_helices" }):
        '''
        VEP anno string is delimited by '&', each domain is formatted generally as 'db:id',
        '''
        if not isinstance(vep_anno_str, str):
            self.logger.warning(f"VEP anno string is not a string: {vep_anno_str}")
            return None
        vep_domains = vep_anno_str.split('&')
        norm_domains = set([])
        for vep_domain in vep_domains:
            vep_anno_pack = vep_domain.split(':', 1)
            if len(vep_anno_pack) < 2:
                self.logger.warning(f"VEP anno domain {vep_domain} is not formatted in a standard format: {vep_anno_pack}")
                continue
            db_name, db_id = vep_anno_pack
            normalized_pack = self.normalize_vep_to_interpro(db_name, db_id)
            if normalized_pack is None:
                self.logger.warning(f"No InterPro entries found for VEP anno domain {db_name}:{db_id}")
                continue
            norm_domains.add(normalized_pack)
        if len(norm_domains) == 0:
            self.logger.warning(f"No InterPro entries found for VEP anno {vep_anno_str}")
            vep_domains = [vep_domain for vep_domain in vep_domains if vep_domain.split(':', 1)[0] in vep_pred_domains]
            return {"vep_domains": vep_domains, "interpro_entries": None}
        else:
            self.logger.debug(f"Found {len(norm_domains)} InterPro XML normalized domain names for VEP anno {vep_anno_str}: {norm_domains}")
            '''
            InterProMapValue = (
                            ipr_id,
                            ipr_type,
                            ipr_short_name,
                            ipr_name,
                            go_terms_list,
                            abstract_text
                        )
            where go_terms_list is a list of tuples (go_id, go_category, go_description)
            '''
            interpro_entry_list = []
            for norm_pack in norm_domains:
                interpro_entry = interpro_entry_map_dict.get(norm_pack)
                if interpro_entry is None:
                    self.logger.warning(f"No InterPro entry found for normalized domain {norm_pack}")
                    continue
                else:
                    self.logger.debug(f"Found {len(interpro_entry)} InterPro entries for normalized domain {norm_pack}: {interpro_entry}")
                    interpro_entry_list.extend(interpro_entry)
            vep_domains = [vep_domain for vep_domain in vep_domains if vep_domain.split(':', 1)[0] in vep_pred_domains]
            return {"vep_domains": vep_domains, "interpro_entries": interpro_entry_list}
    

    def load_curated_functional_domains(self, tsv_file_path: str = os.path.join(os.path.dirname(self_directory), 'data', 'InterPro', 'curated_InterPro_func_domains.tsv.gz')) -> Set[str]:
        '''
        Loads the curated functional domains from a TSV file.
        '''
        interpro_ipr_func_dict = {}
        try:
            with gzip.open(tsv_file_path, 'rt', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    ipr_id = row[0] # First column is the InterPro ID
                    ipr_type = row[8] # Ninth column is the InterPro functional status (functional or not)
                    interpro_ipr_func_dict[ipr_id] = ipr_type
        except Exception as e:
            self.logger.exception(f"Error loading curated functional domains from {tsv_file_path}: {e}")
            raise
        self._interpro_ipr_func_dict = interpro_ipr_func_dict
        
    
    def interpret_functionality(self, vep_anno_str:str, interpro_entry_map_dict: FinalMapping, vep_pred_domains= { "Cleavage_site_(Signalp)",
                                                                                                                    "Coiled-coils_(Ncoils)",
                                                                                                                    "Low_complexity_(Seg)",
                                                                                                                    "Transmembrane_helices" }):
        '''
        domain_info is a dictionary of the form {
            "vep_domains": [vep_domain1, vep_domain2, ...],
            "interpro_entries": [interpro_entry1, interpro_entry2, ...]
        }
        '''
        domain_info = self.query_interpro_entry_vep_anno(vep_anno_str, interpro_entry_map_dict, vep_pred_domains)
        if domain_info is None:
            self.logger.warning(f"No domain info found for VEP anno {vep_anno_str}")
            return np.nan

        vep_domains = domain_info["vep_domains"]
        func_vep_pred_domains = vep_pred_domains - { "Low_complexity_(Seg)" }
        func_vep_domains = [vep_domain for vep_domain in vep_domains if vep_domain.split(':', 1)[0] in func_vep_pred_domains]
        if len(func_vep_domains) > 0:
            self.logger.debug(f"Found {len(func_vep_domains)} functional VEP domains for the query domains {domain_info}")
            return "Functional"

        # Now we need to check the InterPro entries
        interpro_entries = domain_info["interpro_entries"]
        if interpro_entries:
            self.logger.debug(f"Found {len(interpro_entries)} InterPro entries for the query domains {domain_info}")
            if any(DomainNormalizer.is_functional_domain(interpro_entry, self._interpro_ipr_func_dict) for interpro_entry in interpro_entries):
                return "Functional"
            else:
                return "Unknown"
        else:
            self.logger.warning(f"No InterPro entries found for the query domains {domain_info}")
            return np.nan


    @staticmethod
    def is_functional_domain(entry_data: List[Any], interpro_ipr_func_dict: Dict[str, str] = None) -> bool:
        """
        Identifies functional InterPro domains/regions with high accuracy by analyzing
        multiple data fields including type, name, and GO terms. Handles cases
        where Families or Superfamilies describe functional domains.

        Args:
            entry_data: A list representing a single InterPro entry's data:
                        [
                            ipr_id (str),          # Index 0
                            ipr_type (str),        # Index 1
                            ipr_short_name (str),  # Index 2
                            ipr_name (str),        # Index 3
                            go_terms_list (list),  # Index 4
                            abstract_text (str)    # Index 5
                        ]
            interpro_ipr_func_dict: A dictionary mapping InterPro IDs to functional status (curated by Gemini 2.0 flash thinking).
        Returns:
            bool: True if the entry is deemed functional based on evidence, False otherwise.
        """
        # --- Basic Checks ---
        if not isinstance(entry_data, tuple) or len(entry_data) < 5:
            raise ValueError(f"Invalid entry data format: {entry_data}")

        try:
            entry_id = entry_data[0]
            entry_type = entry_data[1]
            short_name = entry_data[2]
            full_name = entry_data[3]
            go_terms = entry_data[4]
        except IndexError:
            raise ValueError(f"Entry data missing expected elements: {entry_data}")

        if interpro_ipr_func_dict is not None:
            curated_result = interpro_ipr_func_dict.get(entry_id, None)
            if curated_result is not None:
                return curated_result == "yes"

        # --- RULE 1: GO TERM-BASED (Strongest Evidence) ---
        # Check GO terms first, as they provide direct functional annotation.
        if isinstance(go_terms, list) and go_terms:
            for go_term in go_terms:
                if isinstance(go_term, list) and len(go_term) >= 3:
                    go_category = go_term[1]
                    go_description = go_term[2]

                    # Rule 1.1: Any molecular function GO term is strong evidence
                    if go_category == "molecular_function":
                        # logger.debug(f"Functional GO (molecular_function) found for {entry_id}")
                        return True

                    # Rule 1.2: Specific biological processes implying function
                    if go_category == "biological_process":
                        functional_processes = {
                            'catalysis', 'transport', 'signaling', 'regulation',
                            'biosynthetic', 'metabolic', 'synthesis', 'activity',
                            'binding', 'interaction', 'transcription', 'replication',
                            'repair', 'modification' # Added more specific process types
                        }
                        desc_lower = go_description.lower()
                        if any(proc in desc_lower for proc in functional_processes):
                            # logger.debug(f"Functional GO (biological_process: {go_description}) found for {entry_id}")
                            return True

        # --- RULE 2: NAME-BASED (Strong Supporting Evidence) ---
        # Check names if GO terms didn't provide a clear answer.
        # Crucial for entries like Families/Superfamilies describing domains.
        functional_keywords = {
            # Activities/Roles
            'enzyme', 'activity', 'binding', 'receptor', 'kinase', 'transferase',
            'hydrolase', 'isomerase', 'ligase', 'lyase', 'reductase', 'oxidase',
            'synthase', 'synthetase', 'phosphatase', 'channel', 'transporter',
            'carrier', 'permease', 'catalytic', 'regulatory', 'effector',
            'interaction', 'recognition', 'transcription', 'replication', 'repair',
            'modification', 'inhibitor', 'activator',
            # Specific functional domain types mentioned in names
            'homeobox', 'helix-loop-helix', 'zinc finger', 'zn finger', 'sh2', 'sh3',
            'wd40', 'ankyrin', 'bzip', 'bromodomain', 'chromodomain',
            # General terms often indicating function
            'active site', 'binding site', 'catalytic site', 'substrate', 'ligand'
        }
        # Include common suffixes/patterns indicating domains or function
        functional_patterns = [
            r'ase$', r'_dom$', r'domain$', r'motif$', r'finger$', r'box$', # Common functional suffixes/types
            r'active_site', r'binding_site', r'catalytic_site',
            r'dna-binding', r'rna-binding', r'protein-binding' # Common binding descriptions
        ]

        name_text = (short_name + " " + full_name).lower()

        if any(keyword in name_text for keyword in functional_keywords):
            # logger.debug(f"Functional keyword found in name for {entry_id}")
            return True
        if any(re.search(pattern, name_text) for pattern in functional_patterns):
            # logger.debug(f"Functional pattern found in name for {entry_id}")
            return True

        # --- RULE 3: TYPE-BASED (Weakest Evidence / Fallback) ---
        # Only rely on type if GO and Name checks didn't find anything conclusive.
        # These types are inherently functional or represent specific functional regions.
        if entry_type in ["Active_site", "Binding_site", "Conserved_site", "PTM"]: # Post-translational modification site
             # logger.debug(f"Functional type ({entry_type}) found for {entry_id}")
             return True

        # If it's a "Domain" type but didn't trigger GO or Name rules, we might still
        # consider it functional depending on desired sensitivity vs accuracy.
        # For higher accuracy, we require positive evidence from GO/Name.
        # if entry_type == "Domain":
        #     logger.debug(f"Type 'Domain' found for {entry_id}, but no strong GO/Name evidence. Classified as non-functional for accuracy.")
        #     # return True # Uncomment for higher sensitivity

        # --- Default: Not Functional ---
        # If none of the above rules triggered, assume non-functional based on available evidence.
        # logger.debug(f"No strong functional evidence found for {entry_id} (Type: {entry_type}). Classified as non-functional.")
        return False

# --- Standalone Function to Parse InterPro XML for the Full Mapping ---
# (This function remains the same as in the previous response, including
# the removal of unused variables and correction of return type hint)
def parse_interpro_xml_for_mapping(xml_file_path: str) -> FinalMapping:
    """
    Parses a gzipped InterPro XML file (interpro.xml.gz) to extract mappings
    from member database signatures to InterPro entries, including GO terms
    and abstract text. (Single-process version).
    """
    logger.info(f"Starting single-process parsing of {xml_file_path} for FULL mapping...")
    start_time = time.time()
    member_to_interpro_map: FinalMapping = defaultdict(list)
    interpro_entry_count = 0

    try:
        with gzip.open(xml_file_path, 'rb') as f:
            context = ET.iterparse(f, events=('end',))
            _, root = next(context)

            for event, elem in context:
                if elem.tag == 'interpro':
                    interpro_entry_count += 1
                    ipr_id = None
                    try:
                        ipr_id = elem.get('id')
                        ipr_type = elem.get('type')
                        ipr_short_name = elem.get('short_name')
                        name_element = elem.find('name')
                        ipr_name = name_element.text.strip() if name_element is not None and name_element.text else ""

                        go_terms_list: List[GOTermInfo] = []
                        class_list = elem.find('class_list')
                        if class_list is not None:
                            for classification in class_list.findall('classification'):
                                if classification.get('class_type') == 'GO':
                                    go_id = classification.get('id')
                                    cat_elem = classification.find('category')
                                    desc_elem = classification.find('description')
                                    go_category = cat_elem.text.strip() if cat_elem is not None and cat_elem.text else ""
                                    go_description = desc_elem.text.strip() if desc_elem is not None and desc_elem.text else ""
                                    if go_id:
                                        go_terms_list.append((go_id, go_category, go_description))

                        abstract_text = ""
                        abstract_element = elem.find('abstract')
                        if abstract_element is not None:
                            abstract_text = " ".join(t for t in abstract_element.itertext() if t).strip()
                            abstract_text = ' '.join(abstract_text.split())

                        map_value: InterProMapValue = (
                            ipr_id or "UnknownIPR",
                            ipr_type or "UnknownType",
                            ipr_short_name or "UnknownShortName",
                            ipr_name,
                            go_terms_list,
                            abstract_text
                        )

                        member_list = elem.find('member_list')
                        if member_list is not None and ipr_id is not None:
                            for db_xref in member_list.findall('db_xref'):
                                db_name = db_xref.get('db')
                                db_key = db_xref.get('dbkey')
                                if db_name and db_key:
                                    map_key = (db_name.upper(), db_key)
                                    member_to_interpro_map[map_key].append(map_value)

                    except Exception as inner_e:
                         logger.exception(f"Error processing details for InterPro entry near count {interpro_entry_count} (ID: {ipr_id})")
                    finally:
                         elem.clear()

                    if interpro_entry_count % 10000 == 0:
                        elapsed = time.time() - start_time
                        logger.info(f"Processed {interpro_entry_count} entries (full map)... ({elapsed:.1f}s)")

        end_time = time.time()
        logger.info(f"\nFinished FULL mapping parsing.")
        logger.info(f"Processed {interpro_entry_count} InterPro entries.")
        logger.info(f"Created mapping dictionary with {len(member_to_interpro_map)} unique member signature keys.")
        logger.info(f"Total time: {end_time - start_time:.2f} seconds.")
        return dict(member_to_interpro_map)

    except FileNotFoundError:
        logger.error(f"Error: File not found at {xml_file_path}")
        return {}
    except ET.ParseError as e:
        logger.exception(f"Error: Failed to parse XML file")
        return {}
    except Exception as e:
        logger.exception(f"An unexpected error occurred during file processing")
        return {}


# --- Main Execution Block (Sequential Workflow) ---
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Collect VEP/InterPro domain formats, generate comparison report, AND create full InterPro mapping pickle file in one run."
    )
    parser.add_argument('--vcf', required=True, help='Path to VCF file with VEP annotations')
    parser.add_argument('--xml', required=True, help='Path to interpro.xml.gz file')
    parser.add_argument('--report', required=True, help='Output path for the comparison report text file')
    parser.add_argument('--mapping_output', required=True, help='Output pickle file for the full InterPro mapping dictionary')

    args = parser.parse_args()

    # --- Step 1: Data Collection and Report Generation using DomainNormalizer ---
    logger.info("--- Running Step 1: Data Collection & Report Generation ---")
    normalizer = DomainNormalizer() # Instantiate class

    logger.info("--- Collecting VEP Domain Formats ---")
    try:
        normalizer.parse_vcf_file(args.vcf)
    except Exception as e:
        logger.error(f"Stopping: Error during VCF processing - {e}", exc_info=True)
        sys.exit(1)

    logger.info("\n--- Collecting InterPro Domain Formats (First Pass) ---")
    try:
        # This only collects names and examples for the report
        normalizer.parse_interpro_xml_formats(args.xml)
    except Exception as e:
        logger.error(f"Stopping: Error during InterPro XML format collection - {e}", exc_info=True)
        sys.exit(1)

    logger.info("\n--- Generating Comparison Report ---")
    try:
        # This generates the report based on data collected above
        normalizer.generate_comparison_report(args.report)
        logger.info(f"\n*** Comparison report generated successfully. Please review and curate mappings. ***")
        logger.info(f"Report saved to: {args.report}")
        logger.info("Implement curated mappings in DomainNormalizer class before using normalization methods.")
    except Exception as e:
        logger.error(f"Stopping: Failed to generate comparison report - {e}", exc_info=True)
        sys.exit(1)

    # --- Step 2: Create Full Mapping Dictionary using Standalone Parser ---
    logger.info("\n--- Running Step 2: Create Full InterPro Mapping Pickle ---")
    logger.warning("Parsing InterPro XML file again to build the full mapping dictionary (including GO/Abstracts). This is inefficient but follows the current script structure.")
    try:
        # This parses the XML again to build the dictionary needed for analysis
        mapping_data = parse_interpro_xml_for_mapping(args.xml)
    except Exception as e:
        logger.error(f"Stopping: Error during full InterPro XML parsing for mapping - {e}", exc_info=True)
        sys.exit(1)

    # --- Step 3: Save the Full Mapping Data ---
    if mapping_data:
        logger.info(f"\n--- Step 3: Saving Full Mapping Data ---")
        logger.info(f"Saving mapping data ({len(mapping_data)} keys) to {args.mapping_output}...")
        try:
            with open(args.mapping_output, 'wb') as f_out:
                pickle.dump(mapping_data, f_out, protocol=pickle.HIGHEST_PROTOCOL)
            logger.info(f"Save successful to {args.mapping_output}")
        except Exception as e:
            logger.error(f"Error saving mapping data - {e}", exc_info=True)
            # Decide if this error is fatal or not
            # sys.exit(1) # Optionally exit if saving fails
    else:
        logger.warning("No InterPro mapping data was generated in Step 2, likely due to errors.")

    logger.info("\n--- Script Finished ---")