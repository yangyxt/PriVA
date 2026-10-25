#!/usr/bin/env bash
# Adapted workflow ordering from CADD-scripts v1.7.2.
# Copyright (c) University of Washington, Hudson-Alpha Institute for
# Biotechnology and Berlin Institute of Health 2013-2020. All rights reserved.
# Permission is hereby granted, to all non-commercial users and licensees of
# CADD, to use, copy, modify, merge, and distribute this software subject to
# retaining the copyright and permission notices. See CADD-LICENSE for the
# complete upstream notice and disclaimer.
set -euo pipefail
export LC_ALL=C

die() { printf 'PriVA CADD annotate: %s\n' "$*" >&2; exit 1; }
[[ $# == 3 ]] || die 'usage: annotate.sh INPUT.vcf[.gz] OUTPUT.anno.tsv.gz WORKDIR'
input=$1
output=$2
workdir=$3
[[ -f "$input" ]] || die "input is not a regular file: $input"
[[ ! -e "$output" && ! -L "$output" ]] || die "output already exists: $output"
[[ -d "$workdir" ]] || die "work directory does not exist: $workdir"

require_env() { [[ -n "${!1:-}" ]] || die "required environment variable is missing: $1"; }
require_file() { [[ -f "$1" && -r "$1" ]] || die "required file is not readable: $1"; }
require_program() { command -v -- "$1" >/dev/null 2>&1 || die "executable is unavailable: $1"; }

CADD_ANNOTATION_INPUT=${CADD_ANNOTATION_INPUT:-shared-vep}
case "$CADD_ANNOTATION_INPUT" in shared-vep|reference-vcf) ;; *) die "unsupported annotation input mode" ;; esac
if [[ "$CADD_ANNOTATION_INPUT" == reference-vcf ]]; then
    for key in CADD_VEP CADD_VEP_CACHE CADD_ENSEMBL_DB; do require_env "$key"; done
    [[ "$CADD_ENSEMBL_DB" =~ ^[1-9][0-9]*$ ]] || die "CADD_ENSEMBL_DB must be a positive release number"
    require_program "$CADD_VEP"
    [[ -d "$CADD_VEP_CACHE" ]] || die "VEP cache directory is missing: $CADD_VEP_CACHE"
fi

for key in CADD_SOURCE_ROOT CADD_PYTHON2 CADD_ASSEMBLY CADD_ESM_PYTHON CADD_ESM_MODEL_DIR CADD_ESM_MODELS \
    CADD_ESM_TRANSCRIPTS CADD_REGSEQ_PYTHON CADD_REGSEQ_MODEL \
    CADD_REGSEQ_WEIGHTS CADD_REGSEQ_REFERENCE CADD_REGSEQ_GENOME \
    CADD_REFERENCE_CONFIG; do
    require_env "$key"
done
case "$CADD_ASSEMBLY" in
    GRCh37) ;;
    GRCh38)
        for key in CADD_MMSPLICE_PYTHON CADD_MMSPLICE_GTF CADD_REFERENCE_FASTA; do
            require_env "$key"
        done
        ;;
    *) die 'CADD_ASSEMBLY must be GRCh37 or GRCh38' ;;
esac

CADD_BGZIP=${CADD_BGZIP:-bgzip}
CADD_TABIX=${CADD_TABIX:-tabix}
CADD_ESM_BATCH_SIZE=${CADD_ESM_BATCH_SIZE:-1}
[[ "$CADD_ESM_BATCH_SIZE" =~ ^[1-9][0-9]*$ ]] || die 'CADD_ESM_BATCH_SIZE must be positive'
for program in git gzip sort uniq mktemp cat "$CADD_PYTHON2" \
    "$CADD_ESM_PYTHON" "$CADD_REGSEQ_PYTHON" "$CADD_BGZIP"; do
    require_program "$program"
done

expected_commit=1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866
actual_commit=$(git -C "$CADD_SOURCE_ROOT" rev-parse HEAD) || die 'upstream source must be a Git checkout'
[[ "$actual_commit" == "$expected_commit" ]] || die "upstream source revision differs from $expected_commit"
git -C "$CADD_SOURCE_ROOT" diff --quiet HEAD -- src config || die 'upstream tracked src/config differs from the pinned revision'
scripts="$CADD_SOURCE_ROOT/src/scripts"
for file in "$scripts/VCF2vepVCF.py" "$scripts/annotateVEPvcf.py" \
    "$scripts/lib/tools/esmScore/esmScore_missense_av_fast.py" \
    "$scripts/lib/tools/esmScore/esmScore_frameshift_av.py" \
    "$scripts/lib/tools/esmScore/esmScore_inFrame_av.py" \
    "$scripts/lib/tools/regulatorySequence/predictVariants.py" \
    "$CADD_ESM_TRANSCRIPTS" "$CADD_REGSEQ_MODEL" "$CADD_REGSEQ_WEIGHTS" \
    "$CADD_REGSEQ_REFERENCE" "$CADD_REGSEQ_REFERENCE.fai" \
    "$CADD_REGSEQ_GENOME" "$CADD_REFERENCE_CONFIG"; do
    require_file "$file"
done
[[ -d "$CADD_ESM_MODEL_DIR" ]] || die "ESM hub directory is missing: $CADD_ESM_MODEL_DIR"

# torch.hub.set_dir expects the parent of checkpoints/. Check local model
# files before calling the upstream loader, which otherwise may download them.
[[ "$CADD_ESM_MODELS" != ,* && "$CADD_ESM_MODELS" != *, && "$CADD_ESM_MODELS" != *,,* ]] \
    || die 'CADD_ESM_MODELS must be a comma-separated list without empty items'
IFS=',' read -r -a model_names <<< "$CADD_ESM_MODELS"
model_args=()
for model in "${model_names[@]}"; do
    [[ "$model" =~ ^esm1v_t33_650M_UR90S_[1-5]$ ]] || die "unsupported ESM-1v model name: $model"
    require_file "$CADD_ESM_MODEL_DIR/checkpoints/$model.pt"
    model_args+=(--model "$model")
done

if [[ "$CADD_ASSEMBLY" == GRCh38 ]]; then
    for program in "$CADD_MMSPLICE_PYTHON" "$CADD_TABIX" grep; do require_program "$program"; done
    for file in "$scripts/lib/tools/MMSplice.py" "$CADD_MMSPLICE_GTF" \
        "$CADD_REFERENCE_FASTA" "$CADD_REFERENCE_FASTA.fai"; do require_file "$file"; done
fi

stage=$(mktemp -d "$workdir/annotation.XXXXXXXX")
# Keep all intermediate data and tool temporary files in the caller's scratch.
export TMPDIR="$stage"
export CADD="$CADD_SOURCE_ROOT"
export PYTHONDONTWRITEBYTECODE=1
read_input() {
    case "$input" in *.gz) gzip -cd -- "$input" ;; *) cat -- "$input" ;; esac
}

# Production consumes the qualified CADD view from PriVA's single official VEP
# run. Re-running VEP is available only as explicit reference validation.
if [[ "$CADD_ANNOTATION_INPUT" == reference-vcf ]]; then
    read_input | "$CADD_PYTHON2" "$scripts/VCF2vepVCF.py" \
        | sort -k1,1 -k2,2n -k4,4 -k5,5 | uniq > "$stage/prepared.vcf"
    "$CADD_VEP" --quiet --cache --offline --dir "$CADD_VEP_CACHE" \
        --buffer 1000 --no_stats --species homo_sapiens \
        --db_version="$CADD_ENSEMBL_DB" --assembly "$CADD_ASSEMBLY" \
        --format vcf --regulatory --sift b --polyphen b --per_gene --ccds --domains \
        --numbers --canonical --total_length --vcf --force_overwrite --output_file STDOUT \
        < "$stage/prepared.vcf" | "$CADD_BGZIP" -c > "$stage/vep.vcf.gz"
else
    read_input | "$CADD_BGZIP" -c > "$stage/vep.vcf.gz"
fi

"$CADD_ESM_PYTHON" "$scripts/lib/tools/esmScore/esmScore_missense_av_fast.py" \
    --input "$stage/vep.vcf.gz" --transcripts "$CADD_ESM_TRANSCRIPTS" \
    --model-directory "$CADD_ESM_MODEL_DIR" "${model_args[@]}" \
    --output "$stage/esm_missense.vcf.gz" --batch-size "$CADD_ESM_BATCH_SIZE" \
    > "$stage/esm_missense.log" 2>&1
"$CADD_ESM_PYTHON" "$scripts/lib/tools/esmScore/esmScore_frameshift_av.py" \
    --input "$stage/esm_missense.vcf.gz" --transcripts "$CADD_ESM_TRANSCRIPTS" \
    --model-directory "$CADD_ESM_MODEL_DIR" "${model_args[@]}" \
    --output "$stage/esm_frameshift.vcf.gz" --batch-size "$CADD_ESM_BATCH_SIZE" \
    > "$stage/esm_frameshift.log" 2>&1
"$CADD_ESM_PYTHON" "$scripts/lib/tools/esmScore/esmScore_inFrame_av.py" \
    --input "$stage/esm_frameshift.vcf.gz" --transcripts "$CADD_ESM_TRANSCRIPTS" \
    --model-directory "$CADD_ESM_MODEL_DIR" "${model_args[@]}" \
    --output "$stage/esm.vcf.gz" --batch-size "$CADD_ESM_BATCH_SIZE" \
    > "$stage/esm_inframe.log" 2>&1
"$CADD_REGSEQ_PYTHON" "$scripts/lib/tools/regulatorySequence/predictVariants.py" \
    --variants "$stage/esm.vcf.gz" --model "$CADD_REGSEQ_MODEL" \
    --weights "$CADD_REGSEQ_WEIGHTS" --reference "$CADD_REGSEQ_REFERENCE" \
    --genome "$CADD_REGSEQ_GENOME" --output "$stage/regseq.vcf.gz" \
    > "$stage/regseq.log" 2>&1

annotation_vcf="$stage/regseq.vcf.gz"
if [[ "$CADD_ASSEMBLY" == GRCh38 ]]; then
    "$CADD_TABIX" -p vcf "$stage/regseq.vcf.gz"
    KERAS_BACKEND=tensorflow "$CADD_MMSPLICE_PYTHON" "$scripts/lib/tools/MMSplice.py" \
        -i "$stage/regseq.vcf.gz" -g "$CADD_MMSPLICE_GTF" -f "$CADD_REFERENCE_FASTA" \
        | grep -v '^Variant(CHROM=' | "$CADD_BGZIP" -c > "$stage/mmsplice.vcf.gz"
    annotation_vcf="$stage/mmsplice.vcf.gz"
fi
gzip -cd -- "$annotation_vcf" \
    | "$CADD_PYTHON2" "$scripts/annotateVEPvcf.py" -c "$CADD_REFERENCE_CONFIG" \
    | gzip -c > "$stage/annotations.tsv.gz"
gzip -t -- "$stage/annotations.tsv.gz"
# The task adapter validates and publishes this staged output.
(set -o noclobber; cat -- "$stage/annotations.tsv.gz" > "$output")
