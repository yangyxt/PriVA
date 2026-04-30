#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

usage() {
    echo "Usage: $0 <hg19|hg38|GRCh37|GRCh38> [config_path] [cores]" >&2
    exit 1
}

[[ $# -ge 1 ]] || usage

assembly_arg="${1}"
config_arg="${2:-}"
cores="${3:-${PRIVA_SMOKE_CORES:-2}}"

assembly_norm=""
case "${assembly_arg}" in
    hg19|GRCh37) assembly_norm="hg19" ;;
    hg38|GRCh38) assembly_norm="hg38" ;;
    *) echo "ERROR: Unsupported assembly: ${assembly_arg}. Use hg19/hg38/GRCh37/GRCh38" >&2; exit 1 ;;
esac

if [[ -z "${config_arg}" ]]; then
    config_arg="example/configs/smoke_demo20_${assembly_norm}.yaml"
fi

to_abs_path() {
    local p="$1"
    if [[ "${p}" = /* ]]; then
        printf '%s\n' "${p}"
    else
        printf '%s\n' "${REPO_ROOT}/${p}"
    fi
}

config_file="$(to_abs_path "${config_arg}")"
[[ -f "${config_file}" ]] || { echo "ERROR: Config file not found: ${config_file}" >&2; exit 1; }
[[ "${cores}" =~ ^[0-9]+$ ]] && [[ "${cores}" -ge 1 ]] || { echo "ERROR: Invalid cores: ${cores}" >&2; exit 1; }

for tool in bash yq bcftools tabix samtools snakemake python; do
    command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: Required tool not found in PATH: ${tool}" >&2; exit 1; }
done

get_cfg() {
    local key="$1"
    local v
    if yq e ".${key}" "${config_file}" >/dev/null 2>&1; then
        v="$(yq e ".${key}" "${config_file}")"
    else
        v="$(yq -r ".${key}" "${config_file}")"
    fi
    [[ "${v}" != "null" ]] && [[ -n "${v}" ]] || { echo "ERROR: Missing required config key: ${key}" >&2; exit 1; }
    printf '%s\n' "${v}"
}

check_file() {
    local p="$1"
    [[ -f "${p}" ]] || { echo "ERROR: Missing file: ${p}" >&2; exit 1; }
}

check_dir() {
    local p="$1"
    [[ -d "${p}" ]] || { echo "ERROR: Missing directory: ${p}" >&2; exit 1; }
}

check_vcf_index() {
    local vcf="$1"
    [[ -f "${vcf}.tbi" || -f "${vcf}.csi" ]] || { echo "ERROR: Missing VCF index (.tbi/.csi): ${vcf}" >&2; exit 1; }
}

validate_csq_arity() {
    local vcf="$1"
    local label="$2"
    python - "${vcf}" "${label}" <<'PY'
import gzip
import re
import sys

vcf_path = sys.argv[1]
label = sys.argv[2]

opener = gzip.open if vcf_path.endswith(".gz") else open
csq_field_count = None
mismatch_count = 0
csq_value_count = 0
bad_examples = []

with opener(vcf_path, "rt") as fh:
    for line in fh:
        if line.startswith("##INFO=<ID=CSQ"):
            m = re.search(r"Format: ([^\">]+)", line)
            if m:
                csq_field_count = len(m.group(1).split("|"))
        elif line.startswith("#"):
            continue
        else:
            parts = line.rstrip("\n").split("\t")
            info = parts[7]
            for kv in info.split(";"):
                if not kv.startswith("CSQ="):
                    continue
                for entry in kv[4:].split(","):
                    if entry in ("", "."):
                        continue
                    csq_value_count += 1
                    field_count = entry.count("|") + 1
                    if csq_field_count is not None and field_count != csq_field_count:
                        mismatch_count += 1
                        if len(bad_examples) < 5:
                            bad_examples.append(
                                f"{parts[0]}:{parts[1]} {parts[3]}>{parts[4]} value_fields={field_count} expected={csq_field_count} CSQ={entry[:160]}"
                            )

if csq_field_count is None:
    print(f"ERROR: {label} missing INFO/CSQ header: {vcf_path}", file=sys.stderr)
    sys.exit(1)

if mismatch_count > 0:
    print(
        f"ERROR: {label} has malformed INFO/CSQ values in {vcf_path}. "
        f"Header declares {csq_field_count} fields but found {mismatch_count} mismatched CSQ entries "
        f"out of {csq_value_count}.",
        file=sys.stderr,
    )
    for ex in bad_examples:
        print(f"ERROR: Example bad CSQ: {ex}", file=sys.stderr)
    sys.exit(1)

print(f"Validated CSQ arity for {label}: {vcf_path} (header_fields={csq_field_count}, csq_values={csq_value_count})")
PY
}

config_assembly="$(get_cfg assembly)"
if [[ "${config_assembly}" != "${assembly_norm}" ]] && \
   ! { [[ "${assembly_norm}" == "hg19" ]] && [[ "${config_assembly}" == "GRCh37" ]]; } && \
   ! { [[ "${assembly_norm}" == "hg38" ]] && [[ "${config_assembly}" == "GRCh38" ]]; }; then
    echo "ERROR: Assembly mismatch. Requested=${assembly_norm}, config assembly=${config_assembly}" >&2
    exit 1
fi

input_vcf="$(to_abs_path "$(get_cfg input_vcf)")"
output_dir="$(to_abs_path "$(get_cfg output_dir)")"
hub_vcf_file="$(to_abs_path "$(get_cfg hub_vcf_file)")"
hub_cadd_file="$(to_abs_path "$(get_cfg hub_cadd_file)")"
ref_genome="$(to_abs_path "$(get_cfg ref_genome)")"
gnomad_vcf_chrX="$(to_abs_path "$(get_cfg gnomad_vcf_chrX)")"
clinvar_vcf="$(to_abs_path "$(get_cfg clinvar_vcf)")"
control_vcf="$(to_abs_path "$(get_cfg control_vcf)")"
clinvar_aa_stat="$(to_abs_path "$(get_cfg clinvar_aa_stat)")"
clinvar_splice_stat="$(to_abs_path "$(get_cfg clinvar_splice_stat)")"
clinvar_patho_af_stat="$(to_abs_path "$(get_cfg clinvar_patho_af_stat)")"
clinvar_patho_exon_af_stat="$(to_abs_path "$(get_cfg clinvar_patho_exon_af_stat)")"
alphamissense_tranx_domain_map="$(to_abs_path "$(get_cfg alphamissense_tranx_domain_map)")"
alphamissense_intolerant_domains="$(to_abs_path "$(get_cfg alphamissense_intolerant_domains)")"
all_intolerant_domains="$(to_abs_path "$(get_cfg all_intolerant_domains)")"
pm1_regions_pkl="$(to_abs_path "$(get_cfg pm1_regions_pkl)")"
interpro_mapping_pickle="$(to_abs_path "$(get_cfg interpro_mapping_pickle)")"
clingen_map="$(to_abs_path "$(get_cfg clingen_map)")"
gene_dosage_sensitivity="$(to_abs_path "$(get_cfg gene_dosage_sensitivity)")"
alphamissense_vcf="$(to_abs_path "$(get_cfg alphamissense_vcf)")"
tmp_dir="$(to_abs_path "$(get_cfg tmp_dir)")"

repeat_region_file_name="$(get_cfg repeat_region_file_name)"
repeat_region_file="${REPO_ROOT}/data/repeats/${repeat_region_file_name}"
clinvar_gene_stat_basename="$(basename "$(get_cfg clinvar_gene_stat)")"
clinvar_gene_stat_file="${REPO_ROOT}/data/ClinVar/${clinvar_gene_stat_basename}"

if pgrep -af "snakemake.*${config_file}" >/dev/null 2>&1; then
    echo "ERROR: Existing snakemake process already using ${config_file}" >&2
    exit 1
fi
if [[ -e "${output_dir}/.snakemake/lock" ]]; then
    echo "ERROR: Detected existing Snakemake lock at ${output_dir}/.snakemake/lock" >&2
    exit 1
fi

check_file "${input_vcf}"
check_vcf_index "${input_vcf}"
check_file "${hub_vcf_file}"
check_vcf_index "${hub_vcf_file}"
check_file "${hub_cadd_file}"
check_file "${gnomad_vcf_chrX}"
check_vcf_index "${gnomad_vcf_chrX}"
check_file "${clinvar_vcf}"
check_vcf_index "${clinvar_vcf}"
check_file "${control_vcf}"
check_vcf_index "${control_vcf}"
check_file "${clinvar_aa_stat}"
check_file "${clinvar_splice_stat}"
check_file "${clinvar_patho_af_stat}"
check_file "${clinvar_patho_exon_af_stat}"
check_file "${alphamissense_tranx_domain_map}"
check_file "${alphamissense_intolerant_domains}"
check_file "${all_intolerant_domains}"
check_file "${pm1_regions_pkl}"
check_file "${interpro_mapping_pickle}"
check_file "${clingen_map}"
check_file "${gene_dosage_sensitivity}"
check_file "${alphamissense_vcf}"
check_file "${repeat_region_file}"
check_file "${clinvar_gene_stat_file}"
validate_csq_arity "${input_vcf}" "input_vcf"
validate_csq_arity "${hub_vcf_file}" "hub_vcf_file"

for d in "${output_dir}" "${tmp_dir}"; do
    mkdir -p "${d}"
done

csq_header="$(bcftools view -h "${input_vcf}" | grep '^##INFO=<ID=CSQ' || true)"
[[ -n "${csq_header}" ]] || { echo "ERROR: Input VCF missing INFO/CSQ header: ${input_vcf}" >&2; exit 1; }

required_csq_fields=(
    "SpliceVault_top_events"
    "SpliceAI_pred_DS_AG"
    "SpliceAI_pred_DS_AL"
    "SpliceAI_pred_DS_DG"
    "SpliceAI_pred_DS_DL"
    "am_pathogenicity"
    "am_class"
    "LoF"
    "LoF_filter"
    "LoF_flags"
    "LoF_info"
)
missing_fields=()
for f in "${required_csq_fields[@]}"; do
    [[ "${csq_header}" == *"${f}"* ]] || missing_fields+=("${f}")
done
if [[ ${#missing_fields[@]} -gt 0 ]]; then
    echo "ERROR: CSQ header missing required fields: ${missing_fields[*]}" >&2
    exit 1
fi

missing_hub_keys_count="$(
    comm -23 \
        <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${input_vcf}" | sort -u) \
        <(awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,$2,$3,$4}' "${hub_cadd_file}" | sort -u) \
        | wc -l | tr -d ' '
)"
[[ "${missing_hub_keys_count}" == "0" ]] || {
    echo "ERROR: hub_cadd_file does not fully cover input variants (${missing_hub_keys_count} missing keys): ${hub_cadd_file}" >&2
    exit 1
}

export PRIVA_SMOKE_SKIP_VEP=1
bash "${REPO_ROOT}/scripts/install_utils.sh" reference_genome_install "${config_file}"
check_file "${ref_genome}"
check_file "${ref_genome}.fai"

input_basename="$(basename "${input_vcf}")"
staged_input_vcf="${tmp_dir}/${input_basename}"
cp -f "${input_vcf}" "${staged_input_vcf}"
if [[ -f "${input_vcf}.tbi" ]]; then
    cp -f "${input_vcf}.tbi" "${staged_input_vcf}.tbi"
elif [[ -f "${input_vcf}.csi" ]]; then
    cp -f "${input_vcf}.csi" "${staged_input_vcf}.csi"
fi
check_file "${staged_input_vcf}"
check_vcf_index "${staged_input_vcf}"

staged_base="$(basename "${staged_input_vcf}")"
staged_base="${staged_base%.vcf.gz}"
if [[ "${staged_base}" == *.vcf ]]; then
    staged_base="${staged_base%.vcf}"
fi

collision_candidates=(
    "${output_dir}/${staged_base}.anno.vcf.gz"
    "${output_dir}/${staged_base}.anno.cadd.tsv"
    "${output_dir}/${staged_base}.anno.filtered.vcf.gz"
    "${output_dir}/${staged_base}.anno.filtered.tsv"
    "${output_dir}/${staged_base}.anno.filtered.acmg.tsv"
)
for c in "${collision_candidates[@]}"; do
    if [[ -e "${c}" ]]; then
        echo "ERROR: Output collision detected. Remove or change output_dir before rerun: ${c}" >&2
        exit 1
    fi
done

snakemake \
    --snakefile "${REPO_ROOT}/Snakefile" \
    --configfile "${config_file}" \
    --config "input_vcf=${staged_input_vcf}" "control_vcf=${staged_input_vcf}" \
    --cores "${cores}" \
    --rerun-incomplete

input_base="$(basename "${staged_input_vcf}")"
input_base="${input_base%.vcf.gz}"
if [[ "${input_base}" == *.vcf ]]; then
    input_base="${input_base%.vcf}"
fi

if yq e '.ped_file' "${config_file}" >/dev/null 2>&1; then
    ped_file="$(yq e '.ped_file' "${config_file}")"
else
    ped_file="$(yq -r '.ped_file' "${config_file}")"
fi
if [[ "${ped_file}" != "null" ]] && [[ -n "${ped_file}" ]]; then
    ped_file_abs="$(to_abs_path "${ped_file}")"
else
    ped_file_abs=""
fi

if [[ -n "${ped_file_abs}" ]] && [[ -s "${ped_file_abs}" ]]; then
    mapfile -t fams < <(awk 'NF>0 && $1 !~ /^#/ {print $1}' "${ped_file_abs}" | sort -u)
    [[ ${#fams[@]} -gt 0 ]] || { echo "ERROR: ped_file is non-empty but no families parsed: ${ped_file_abs}" >&2; exit 1; }
    for fam in "${fams[@]}"; do
        check_file "${output_dir}/${input_base}.anno.${fam}.filtered.tsv"
        check_file "${output_dir}/${input_base}.anno.${fam}.filtered.acmg.tsv"
    done
    echo "Smoke test succeeded (${assembly_norm}). Outputs in ${output_dir} for families: ${fams[*]}"
else
    final_tsv="${output_dir}/${input_base}.anno.filtered.tsv"
    final_acmg="${output_dir}/${input_base}.anno.filtered.acmg.tsv"
    check_file "${final_tsv}"
    check_file "${final_acmg}"
    echo "Smoke test succeeded (${assembly_norm}). Outputs: ${final_tsv} and ${final_acmg}"
fi
