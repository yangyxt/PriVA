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

die() { printf 'PriVA CADD calculate: %s\n' "$*" >&2; exit 1; }
[[ $# == 3 ]] || die 'usage: calculate.sh INPUT.anno.tsv[.gz] OUTPUT.scored.tsv WORKDIR'
input=$1
output=$2
workdir=$3
[[ -f "$input" ]] || die "input is not a regular file: $input"
[[ ! -e "$output" && ! -L "$output" ]] || die "output already exists: $output"
[[ -d "$workdir" ]] || die "work directory does not exist: $workdir"
for key in CADD_SOURCE_ROOT CADD_PYTHON2 CADD_IMPUTE_CONFIG CADD_MODEL CADD_PHRED_TABLE; do
    [[ -n "${!key:-}" ]] || die "required environment variable is missing: $key"
done
for program in git gzip mktemp cat "$CADD_PYTHON2"; do
    command -v -- "$program" >/dev/null 2>&1 || die "executable is unavailable: $program"
done
expected_commit=1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866
actual_commit=$(git -C "$CADD_SOURCE_ROOT" rev-parse HEAD) || die 'upstream source must be a Git checkout'
[[ "$actual_commit" == "$expected_commit" ]] || die "upstream source revision differs from $expected_commit"
git -C "$CADD_SOURCE_ROOT" diff --quiet HEAD -- src config || die 'upstream tracked src/config differs from the pinned revision'
scripts="$CADD_SOURCE_ROOT/src/scripts"
for file in "$CADD_IMPUTE_CONFIG" "$CADD_MODEL" "$CADD_PHRED_TABLE" \
    "$scripts/trackTransformation.py" "$scripts/predictSKmodel.py" \
    "$scripts/max_line_hierarchy.py" "$scripts/appendPHREDscore.py"; do
    [[ -f "$file" && -r "$file" ]] || die "required file is not readable: $file"
done

stage=$(mktemp -d "$workdir/calculation.XXXXXXXX")
export TMPDIR="$stage"
export CADD="$CADD_SOURCE_ROOT"
export PYTHONDONTWRITEBYTECODE=1
# Upstream predictSKmodel requires a gzip annotation file for its append input.
# Compression conversion does not alter the annotation bytes or row order.
case "$input" in
    *.gz) gzip -t -- "$input"; annotations="$input" ;;
    *) gzip -c -- "$input" > "$stage/annotations.tsv.gz"; annotations="$stage/annotations.tsv.gz" ;;
esac
gzip -cd -- "$annotations" \
    | "$CADD_PYTHON2" "$scripts/trackTransformation.py" -b \
        -c "$CADD_IMPUTE_CONFIG" -o "$stage/features.csv.gz" --noheader
"$CADD_PYTHON2" "$scripts/predictSKmodel.py" \
    -i "$stage/features.csv.gz" -m "$CADD_MODEL" -a "$annotations" | "$CADD_PYTHON2" "$scripts/max_line_hierarchy.py" --all \
    | "$CADD_PYTHON2" "$scripts/appendPHREDscore.py" -t "$CADD_PHRED_TABLE" \
    > "$stage/scored.tsv"
# All four scientific stages above are original, pinned upstream Python.
(set -o noclobber; cat -- "$stage/scored.tsv" > "$output")
