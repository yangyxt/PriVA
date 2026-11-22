#!/usr/bin/env bash
source "$1/scripts/prioritization_vcf_per_fam.sh" || exit 1
config=$2 step=$3
shift 3
threads=$(read_yaml "$config" threads_per_fam)
case $step in
    combine)
        vcf=$1 cadd=$2 table=${1/.vcf*/.tsv}
        if [[ $(count_vcf_records "$vcf") -lt 1 ]] && check_vcf_validity "$vcf" 0; then : > "$table"; exit 0; fi
        [[ $(wc -l < "$cadd") -le 1 ]] && { echo "The CADD table $cadd is empty" >&2; exit 0; }
        prepare_combined_tab "$vcf" "$cadd" "$threads" "$config" && [[ -f $table ]] ;;
    splicing)
        [[ -s $1 ]] || exit 0
        interpret_splicing_annotations "$1" "$config" "$threads" ;;
    utr)
        [[ -s $1 ]] || exit 0
        interpret_utr_annotations "$1" "$config" "$threads" ;;
    acmg)
        table=$1 family=${2:-}
        [[ -s $table ]] || { : > "${table%.tsv}.acmg.tsv"; exit 0; }
        acmg_assign_table -i "$table" -c "$config" -t "$threads" ${family:+-f "$family"} ;;
    *)
        echo "Unknown prioritization step: $step" >&2; exit 2 ;;
esac
