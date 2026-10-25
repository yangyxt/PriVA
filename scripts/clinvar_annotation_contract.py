#!/usr/bin/env python3
"""Check the annotation fields used by the unchanged ClinVar statistics rules."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import pickle
import hashlib
import os
import sys
import tempfile

import pysam


AA_CSQ_FIELDS = frozenset({
    "Feature", "Feature_type", "Protein_position", "HGVSp", "Consequence",
    "EXON", "INTRON", "HGVSc",
    "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
})
GENE_CSQ_FIELDS = frozenset({
    "SYMBOL", "Gene", "BIOTYPE", "Feature", "DOMAINS", "Consequence",
    "LoF", "NMD", "SpliceVault_top_events", "SpliceVault_SpliceAI_delta",
    "SpliceAI_pred_DS_DL",
})
AF_CSQ_FIELDS = frozenset({"Gene", "Feature", "EXON", "HGVSp", "IMPACT", "Consequence", "SYMBOL"})
AF_INFO_FIELDS = frozenset({"AF_joint", "AF_grpmax_joint"})
# The generic table reader uses these population fields even though the
# gene-statistics rules themselves do not use frequencies.
TABLE_INFO_FIELDS = frozenset({
    "AC_joint", "AN_joint", "AF_joint", "nhomalt_joint", "AC_grpmax_joint",
    "AF_grpmax_joint", "AN_grpmax_joint", "nhomalt_grpmax_joint",
    *(f"{prefix}_{pop}_{sex}"
      for prefix in ("AF_joint",) for pop in
      ("nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining")
      for sex in ("XX", "XY")),
    *(f"nhomalt_joint_{pop}_XX" for pop in
      ("nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining")),
})


def validate_header(path: str | Path, stage: str = "complete") -> dict:
    """Reject incomplete VEP output before a statistics program writes a cache."""
    if stage not in {"aa", "gene", "af", "vep", "complete"}:
        raise ValueError(f"Unknown ClinVar validation stage: {stage}")
    with pysam.VariantFile(str(path)) as vcf:
        info_fields = set(vcf.header.info)
        required_info = {"CSQ", "CLNSIG", "CLNREVSTAT"}
        required_csq = set()
        if stage in {"aa", "vep", "complete"}:
            required_csq.update(AA_CSQ_FIELDS)
        if stage in {"gene", "vep", "complete"}:
            required_csq.update(GENE_CSQ_FIELDS)
            required_info.add("ALLELEID")
        if stage in {"af", "complete"}:
            required_info.update(AF_INFO_FIELDS)
            required_csq.update(AF_CSQ_FIELDS)
        if stage in {"gene", "complete"}:
            required_info.update(TABLE_INFO_FIELDS)
        missing_info = sorted(required_info - info_fields)
        if "CSQ" not in info_fields:
            raise ValueError("ClinVar annotation lacks INFO/CSQ")
        description = vcf.header.info["CSQ"].description
        if "Format: " not in description:
            raise ValueError("ClinVar INFO/CSQ header lacks its Format declaration")
        csq_fields = description.split("Format: ", 1)[1].split("|")
        if len(csq_fields) != len(set(csq_fields)):
            raise ValueError("ClinVar INFO/CSQ header contains duplicate fields")
        missing_csq = sorted(required_csq - set(csq_fields))
        if missing_info or missing_csq:
            raise ValueError(
                f"Incomplete ClinVar {stage} annotation: missing INFO={missing_info}; "
                f"missing CSQ={missing_csq}"
            )
        return {
            "path": str(path), "stage": stage, "success": True,
            "csq_fields": csq_fields, "info_fields": sorted(info_fields),
            "contigs": list(vcf.header.contigs),
        }


def validate_cache(path: str | Path, kind: str) -> dict:
    """Check complete serialization and the existing statistics cache shapes."""
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rb') as handle:
        data = pickle.load(handle)
    if not isinstance(data, dict) or not data:
        raise ValueError(f"ClinVar {kind} cache must contain a nonempty dictionary: {path}")
    if kind == 'aa':
        valid = all(isinstance(positions, dict) and positions and
                    all(isinstance(changes, dict) and changes and
                        all(isinstance(entry, dict) and {'CLNSIG', 'CLNREVSTAT'} <= entry.keys()
                            for entry in changes.values()) for changes in positions.values())
                    for positions in data.values())
    elif kind == 'splice':
        valid = all(isinstance(entries, list) and entries and
                    all(isinstance(entry, dict) and {'chrom', 'pos', 'clinvar_sig', 'splice_ai'} <= entry.keys()
                        for entry in entries) for entries in data.values())
    elif kind == 'gene-af':
        # The original maximum-AF rule can retain None when every AF is zero.
        valid = all(value is None or (isinstance(value, dict) and
                    {'chrom', 'pos', 'af', 'acc_patho_af'} <= value.keys()) for value in data.values())
    elif kind == 'exon-af':
        valid = all(isinstance(exons, dict) and exons and
                    all(isinstance(value, tuple) and len(value) == 5 for value in exons.values())
                    for exons in data.values())
    elif kind == 'gene':
        required = {'gene_symbol', 'ensembl_id', 'clinvar_variants', 'clinvar_pathogenic',
                    'clinvar_not_pathogenic', 'putative_nmd_variants', 'small_aachange',
                    'large_aachange', 'no_aachange', 'transcripts', 'domains'}
        valid = all(isinstance(value, dict) and required <= value.keys() for value in data.values())
    else:
        raise ValueError(f"Unknown ClinVar cache kind: {kind}")
    if not valid:
        raise ValueError(f"ClinVar {kind} cache has an unexpected statistics structure: {path}")
    return {'path': str(path), 'cache_kind': kind, 'success': True, 'top_level_entries': len(data)}


def pair_metadata(path: str | Path) -> dict:
    """Portable metadata identity; never read the VCF or derived payload bytes."""
    path = Path(path).absolute()
    stat = path.stat()
    if not path.is_file() or stat.st_size == 0:
        raise ValueError(f"Missing/non-file/empty ClinVar pair dependency: {path}")
    return dict(path=str(path), target=str(path.resolve()), inode=stat.st_ino,
                size=stat.st_size, mtime_ns=stat.st_mtime_ns, ctime_ns=stat.st_ctime_ns)


def pair_inputs(vcf: str, producer: str) -> dict:
    # pysam prefers CSI when present. Track both existing index candidates so
    # replacing/introducing either cannot silently retain a different fetch view.
    indexes = [str(Path(vcf + suffix)) for suffix in ('.csi', '.tbi')
               if Path(vcf + suffix).exists()]
    if not indexes:
        raise ValueError('ClinVar AA/splice input requires an index')
    return {'vcf': pair_metadata(vcf), 'indexes': [pair_metadata(p) for p in indexes],
            'producer_sha256': hashlib.sha256(Path(producer).read_bytes()).hexdigest(),
            'contract_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def write_pair_json(path: str, value: dict) -> None:
    target = Path(path)
    fd, temporary = tempfile.mkstemp(prefix='.' + target.name + '.', dir=target.parent)
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(value, handle, sort_keys=True)
            handle.write('\n')
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def require_distinct_paths(protected: list[str], *, metadata_paths: list[str] | None = None) -> None:
    """Apply the same read-only path/inode exclusions to ClinVar output pairs."""
    metadata_paths = metadata_paths or []
    resolved = [Path(path).resolve() for path in protected]
    if len(set(resolved)) != len(resolved):
        raise ValueError('ClinVar pair input/output paths must be distinct')
    for path in metadata_paths:
        if Path(path).is_symlink() or Path(path).resolve() in resolved:
            raise ValueError('ClinVar pair metadata must not alias an input or output')
        resolved.append(Path(path).resolve())
    existing = [Path(path) for path in protected + metadata_paths if Path(path).exists()]
    for number, path in enumerate(existing):
        if any(path.samefile(other) for other in existing[number + 1:]):
            raise ValueError('ClinVar pair paths must not share an existing file')


def af_paths_command(argv: list[str]) -> None:
    parser = argparse.ArgumentParser(description='Check ClinVar AF output paths before any writes')
    for name in ('vcf', 'config', 'gene', 'exon'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args(argv)
    # 2026-10-23 CDT: protect the selected VCF, both index names and active
    # configuration before either AF output can create directories or publish.
    require_distinct_paths([args.vcf, args.vcf + '.tbi', args.vcf + '.csi',
                            args.config, args.gene, args.exon])


def pair_command(argv: list[str]) -> None:
    parser = argparse.ArgumentParser(description='ClinVar derived AA/splice completion association')
    parser.add_argument('operation', choices=('pair-capture', 'pair-check', 'pair-complete'))
    for name in ('vcf', 'producer', 'aa', 'splice', 'receipt'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--captured')
    parser.add_argument('--config', help='Protect the caller YAML from output/sidecar aliases')
    args = parser.parse_args(argv)
    protected = [args.vcf, args.vcf + '.tbi', args.vcf + '.csi',
                 args.producer, __file__, args.aa, args.splice]
    if args.config:
        protected.append(args.config)
    metadata_paths = [args.receipt] + ([args.captured] if args.captured else [])
    require_distinct_paths(protected, metadata_paths=metadata_paths)
    inputs = pair_inputs(args.vcf, args.producer)
    if args.operation == 'pair-capture':
        if not args.captured:
            raise ValueError('pair-capture requires --captured')
        write_pair_json(args.captured, inputs)
        return
    state = {'schema': 1, 'identity_policy': 'path-target-inode-size-mtime-ctime; source-content',
             'inputs': inputs, 'aa': pair_metadata(args.aa), 'splice': pair_metadata(args.splice)}
    if args.operation == 'pair-check':
        previous = json.loads(Path(args.receipt).read_text())
        if any(previous.get(key) != value for key, value in state.items()):
            raise ValueError('ClinVar AA/splice completion association is stale')
    else:
        if not args.captured or json.loads(Path(args.captured).read_text()) != inputs:
            raise ValueError('ClinVar AA/splice inputs changed during generation')
        state['producer_path'] = str(Path(args.producer).absolute())  # provenance, not compatibility
        write_pair_json(args.receipt, state)


def main() -> None:
    if len(sys.argv) > 1 and sys.argv[1] == 'af-paths-check':
        try:
            af_paths_command(sys.argv[2:])
        except (OSError, ValueError, TypeError) as error:
            raise SystemExit(str(error))
        return
    if len(sys.argv) > 1 and sys.argv[1].startswith('pair-'):
        try:
            pair_command(sys.argv[1:])
        except (OSError, ValueError, TypeError) as error:
            raise SystemExit(str(error))
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf")
    parser.add_argument("--stage", choices=("aa", "gene", "af", "vep", "complete"),
                        default="complete")
    parser.add_argument("--cache-kind", choices=("aa", "splice", "gene-af", "exon-af", "gene"))
    args = parser.parse_args()
    try:
        result = validate_cache(args.vcf, args.cache_kind) if args.cache_kind else validate_header(args.vcf, args.stage)
    except (OSError, ValueError, TypeError, EOFError, pickle.UnpicklingError) as error:
        parser.exit(1, str(error) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
