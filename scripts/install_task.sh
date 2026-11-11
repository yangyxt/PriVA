#!/usr/bin/env bash
# Runs one install step of conf/install_steps.json as a Nextflow task.
#
#   install_task.sh plan --registry R --config USER --out DIR
#                        [--adopt-config OLD --adopt all|ID,ID] [--run ID,ID]
#   install_task.sh run  --registry R --plan PLAN --step ID --lineage L --folder F
#                        --receipt FILE [--input RECEIPT]... [--stub]
#   install_task.sh adopt --plan PLAN --step ID --lineage L --receipt FILE
#   install_task.sh publish --registry R --plan PLAN --lineage L --folder F
#                           [--input RECEIPT]... [--stub]
#
# plan: writes DIR/base.yaml (your configuration without any key an install
# step writes, so each step resolves its own destinations) and DIR/plan.json
# (per step: code digest, settings, the files you supply, adopted keys).
# run: builds a private configuration whose cache_root is F/data, adds the
# keys of the input receipts, calls the step's installer functions unchanged,
# checks the output keys and writes the receipt. F/data is removed first if a
# failed attempt left it; Nextflow (storeDir) skips a step whose F holds a
# receipt. The finished folder is made read-only unless the step is writable.
# adopt: writes a receipt naming the files of an existing installation.
set -o pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

die() { echo "install_task: $*" >&2; exit 1; }
json_step() { jq -ce --arg id "$2" '.steps[] | select(.id == $id)' "$1" || die "no step $2 in $1"; }
yaml_value() {
    local v
    v=$(yq -r --arg k "$2" '.[$k] // empty | tostring' "$1" 2>/dev/null) || return 1
    [[ $v == null ]] && v=
    printf '%s' "$v"
}
set_key() { # set_key FILE KEY VALUE: replace the key's line or append one
    local file=$1 key=$2 value=$3
    if grep -qE "^${key}[[:space:]]*:" "$file"; then
        local escaped=${value//\\/\\\\}; escaped=${escaped//&/\\&}; escaped=${escaped//|/\\|}
        sed -i -E "s|^(${key}[[:space:]]*:).*|\\1 ${escaped}|" "$file"
    else
        printf '%s: %s\n' "$key" "$value" >> "$file"
    fi
}
identity() { # path size mtime, or "missing"
    local p=$1
    if [[ -e $p ]]; then stat -L -c '%n %s %Y' -- "$p"; else echo "missing $p"; fi
}

# Every function a step's calls reach, by name, then a digest of their text,
# of the script files they name, and of the step's registry entry.
step_digest() {
    local registry=$1 id=$2 entry f w
    entry=$(json_step "$registry" "$id") || return 1
    local -A known=() seen=()
    while read -r _ _ f; do known[$f]=1; done < <(declare -F)
    local -a queue=($(jq -r '.calls[][0]' <<< "$entry"))
    local -a names=()
    while ((${#queue[@]})); do
        f=${queue[0]}; queue=("${queue[@]:1}")
        [[ -n ${seen[$f]:-} ]] && continue
        seen[$f]=1; names+=("$f")
        for w in $(declare -f "$f" | grep -oE '[A-Za-z_][A-Za-z0-9_]*' | sort -u); do
            [[ -n ${known[$w]:-} && -z ${seen[$w]:-} ]] && queue+=("$w")
        done
    done
    {
        printf '%s\n' "$entry"
        for f in $(printf '%s\n' "${names[@]}" | sort); do declare -f "$f"; done
        for w in $(for f in "${names[@]}"; do declare -f "$f"; done | grep -oE '[A-Za-z0-9_.-]+\.(py|sh|R|pl)' | sort -u); do
            [[ -f $SCRIPT_DIR/$w ]] && { printf '%s ' "$w"; sha256sum < "$SCRIPT_DIR/$w"; }
        done
    } | sha256sum | cut -c1-64
}

cmd_plan() {
    local registry config out adopt_config="" adopt="" run=""
    while (($#)); do case $1 in
        --registry) registry=$2; shift 2 ;; --config) config=$2; shift 2 ;; --out) out=$2; shift 2 ;;
        --adopt-config) adopt_config=$2; shift 2 ;; --adopt) adopt=$2; shift 2 ;; --run) run=$2; shift 2 ;;
        *) die "plan: unknown option $1" ;; esac; done
    [[ -r $registry && -r $config && -n $out ]] || die "plan needs --registry, --config and --out"
    mkdir -p "$out"; out=$(cd "$out" && pwd)
    source "$SCRIPT_DIR/install_utils.sh" || die "cannot source install_utils.sh"
    # Keys a step writes, and folder keys, are removed so each step resolves
    # its own destinations under its own cache_root; files you supply stay.
    local base=$out/base.yaml key
    cp "$config" "$base"
    local -A keep=()
    for key in $(jq -r '.steps[].user_inputs[]?' "$registry"); do keep[$key]=1; done
    for key in $(jq -r '(.folder_keys[]), (.steps[] | .outputs[], .optional_outputs[]?)' "$registry" | sort -u); do
        [[ -n ${keep[$key]:-} ]] || sed -i -E "/^${key}[[:space:]]*:/d" "$base"
    done
    local priva_identity
    priva_identity=$(git -C "$BASE_DIR" rev-parse HEAD:rust 2>/dev/null) || \
        priva_identity=$(sha256sum < "${PRIVA_BINARY:?PRIVA_BINARY must name the priva program}" | cut -c1-64)
    local all_keys; all_keys=$(yq -r 'keys[]' "$config")
    local steps="{}" id entry digest settings files adopted pat k v
    for id in $(jq -r '.steps[].id' "$registry"); do
        entry=$(json_step "$registry" "$id")
        digest=$(step_digest "$registry" "$id") || die "digest of $id failed"
        settings="{}"
        for pat in $(jq -r '.settings[]' <<< "$entry"); do
            for k in $all_keys; do
                [[ $k == $pat ]] || continue
                v=$(yaml_value "$config" "$k")
                settings=$(jq -c --arg k "$k" --arg v "$v" '. + {($k): $v}' <<< "$settings")
            done
        done
        files="[]"
        for k in $(jq -r '.user_inputs[]?' <<< "$entry"); do
            v=$(yaml_value "$config" "$k")
            [[ -n $v ]] && files=$(jq -c --arg k "$k" --arg i "$(identity "$v")" '. + [[$k, $i]]' <<< "$files")
        done
        adopted=null
        if [[ -n $adopt_config ]] && { [[ $adopt == all && ",$run," != *",$id,"* ]] || [[ ",$adopt," == *",$id,"* ]]; }; then
            adopted='{"keys":{},"identity":[]}'
            for k in $(jq -r '.outputs[], .optional_outputs[]?' <<< "$entry"); do
                v=$(yaml_value "$adopt_config" "$k")
                [[ -n $v ]] || continue
                adopted=$(jq -c --arg k "$k" --arg v "$v" --arg i "$(identity "$v")" \
                    '.keys += {($k): $v} | .identity += [$i]' <<< "$adopted")
            done
        fi
        steps=$(jq -c --arg id "$id" --arg d "$digest" --argjson s "$settings" --argjson f "$files" \
            --argjson a "$adopted" \
            '. + {($id): {digest: $d, settings: $s, user_files: $f, adopt: $a}}' <<< "$steps")
    done
    jq -n --arg base "$base" --arg pi "$priva_identity" --argjson steps "$steps" \
        '{schema: 1, base_config: $base, priva_identity: $pi, steps: $steps}' > "$out/plan.json"
}

cmd_run() {
    local registry plan id lineage folder receipt stub=0
    local -a inputs=()
    while (($#)); do case $1 in
        --registry) registry=$2; shift 2 ;; --plan) plan=$2; shift 2 ;; --step) id=$2; shift 2 ;;
        --lineage) lineage=$2; shift 2 ;; --folder) folder=$2; shift 2 ;; --receipt) receipt=$2; shift 2 ;;
        --input) inputs+=("$2"); shift 2 ;; --stub) stub=1; shift ;;
        *) die "run: unknown option $1" ;; esac; done
    [[ -r $registry && -r $plan && -n $id && -n $lineage && -n $folder && -n $receipt ]] || die "run needs all options"
    local entry; entry=$(json_step "$registry" "$id") || exit 1
    [[ -e $folder/receipt.json ]] && die "$folder already holds a receipt"
    local work=$folder/data
    if [[ -e $work ]]; then chmod -R u+w "$work" 2>/dev/null; rm -rf -- "$work" || die "cannot clear $work"; fi
    mkdir -p "$work" "$PWD/tmp" || die "cannot create $work"
    source "$SCRIPT_DIR/install_utils.sh" || die "cannot source install_utils.sh"
    local private=$PWD/private.yaml r k v
    cp "$(jq -r .base_config "$plan")" "$private"
    set_key "$private" cache_root "$work"
    set_key "$private" tmp_dir "$PWD/tmp"
    [[ -n ${PRIVA_BINARY:-} ]] && set_key "$private" priva_binary "$PRIVA_BINARY"
    [[ -n ${INSTALL_THREADS:-} ]] && set_key "$private" threads "$INSTALL_THREADS"
    export TMPDIR=$PWD/tmp
    if ((stub == 0)); then
        priva_configure_cache_paths "$private" > resolve.log 2>&1 || { cat resolve.log >&2; die "resolving folders failed"; }
        priva_configure_portable_paths --config "$private" >> resolve.log 2>&1 || { cat resolve.log >&2; die "resolving paths failed"; }
    fi
    for r in "${inputs[@]}"; do
        while IFS=$'\t' read -r k v; do set_key "$private" "$k" "$v"; done < <(jq -r '.keys | to_entries[] | "\(.key)\t\(.value)"' "$r")
    done
    for k in $(jq -r '.copies[]?' <<< "$entry"); do
        v=$(yaml_value "$private" "$k"); [[ -f $v ]] || die "$k ($v) is not a file to copy"
        mkdir -p "$work/copied"; cp -p -- "$v" "$work/copied/"
        local index; for index in "$v.tbi" "$v.csi"; do [[ -f $index ]] && cp -p -- "$index" "$work/copied/"; done
        set_key "$private" "$k" "$work/copied/$(basename "$v")"
    done
    cp "$private" private.before.yaml
    export PRIVA_CACHE_INSTALL_INTERNAL=1 PRIVA_CACHE_INSTALL_BUILD_CONFIG=$private
    local n calls; calls=$(jq -c '.calls[]' <<< "$entry")
    if ((stub)); then
        for k in $(jq -r '.outputs[]' <<< "$entry"); do
            if [[ $k == *_dir || $k == loftee_repo ]]; then mkdir -p "$work/stub/$k"; else mkdir -p "$work/stub"; echo "$id $k" > "$work/stub/$k"; fi
            set_key "$private" "$k" "$work/stub/$k"
        done
    else
        while read -r n; do
            local -a args=()
            while read -r tok; do
                case $tok in
                    '{config}') args+=("$private") ;;
                    '{script_dir}') args+=("$SCRIPT_DIR") ;;
                    '{'*'}') k=${tok#\{}; k=${k%\}}; args+=("$(yaml_value "$private" "$k")") ;;
                    *) args+=("$tok") ;;
                esac
            done < <(jq -r '.[]' <<< "$n")
            echo "== $(date -Is) ${args[*]}"
            "${args[@]}" || die "${args[0]} failed (exit $?)"
        done <<< "$calls"
    fi
    local keys="{}" files="[]" missing=""
    for k in $(jq -r '.outputs[]' <<< "$entry"); do
        v=$(yaml_value "$private" "$k")
        [[ -n $v && -e $v ]] || missing+=" $k"
        keys=$(jq -c --arg k "$k" --arg v "$v" '. + {($k): $v}' <<< "$keys")
    done
    [[ -z $missing ]] || die "outputs missing after the step:$missing"
    for k in $(jq -r '.optional_outputs[]?' <<< "$entry"); do
        v=$(yaml_value "$private" "$k")
        [[ -n $v ]] && keys=$(jq -c --arg k "$k" --arg v "$v" '. + {($k): $v}' <<< "$keys")
    done
    for v in $(jq -r '.[]' <<< "$keys"); do files=$(jq -c --arg i "$(identity "$v")" '. + [$i]' <<< "$files"); done
    local upstream="[]"
    for r in "${inputs[@]}"; do upstream=$(jq -c --slurpfile x "$r" '. + [{id: $x[0].id, lineage: $x[0].lineage}]' <<< "$upstream"); done
    diff private.before.yaml "$private" > config_changes.diff
    [[ $(jq '.writable // false' <<< "$entry") == true ]] || chmod -R a-w "$work"
    jq -n --arg id "$id" --arg name "$(jq -r .name <<< "$entry")" --arg l "$lineage" --argjson keys "$keys" \
        --argjson files "$files" --argjson inputs "$upstream" --arg folder "$folder" --arg host "$(hostname)" \
        --arg done "$(date -Is)" --argjson stub "$stub" \
        '{schema: 1, id: $id, name: $name, lineage: $l, folder: $folder, keys: $keys, files: $files, inputs: $inputs, finished: $done, host: $host, stub: ($stub == 1)}' > "$receipt"
}

# publish: one configuration from your base configuration and every step's
# keys (in registry order, so a later step's key wins), native statistics
# archives converted under F/data, and resources.lock.json naming every
# step's lineage and keys.
cmd_publish() {
    local registry plan lineage folder stub=0
    local -a inputs=()
    while (($#)); do case $1 in
        --registry) registry=$2; shift 2 ;; --plan) plan=$2; shift 2 ;; --lineage) lineage=$2; shift 2 ;;
        --folder) folder=$2; shift 2 ;; --input) inputs+=("$2"); shift 2 ;; --stub) stub=1; shift ;;
        *) die "publish: unknown option $1" ;; esac; done
    [[ -r $registry && -r $plan && -n $lineage && -n $folder ]] || die "publish needs all options"
    local work=$folder/data id r k v
    if [[ -e $work ]]; then chmod -R u+w "$work" 2>/dev/null; rm -rf -- "$work" || die "cannot clear $work"; fi
    mkdir -p "$work" "$PWD/tmp" || die "cannot create $work"
    local merged=$work/config.yaml
    cp "$(jq -r .base_config "$plan")" "$merged"
    set_key "$merged" cache_root "$work"
    local -A by_id=()
    for r in "${inputs[@]}"; do by_id[$(jq -r .id "$r")]=$r; done
    local steps="{}"
    for id in $(jq -r '.steps[].id' "$registry"); do
        r=${by_id[$id]:-}; [[ -n $r ]] || die "no receipt for step $id"
        while IFS=$'\t' read -r k v; do set_key "$merged" "$k" "$v"; done < <(jq -r '.keys | to_entries[] | "\(.key)\t\(.value)"' "$r")
        steps=$(jq -c --slurpfile x "$r" --arg id "$id" '. + {($id): {name: $x[0].name, lineage: $x[0].lineage, adopted: ($x[0].adopted // false), keys: $x[0].keys, receipt: $x[0].folder}}' <<< "$steps")
    done
    [[ -n ${PRIVA_BINARY:-} ]] && set_key "$merged" priva_binary "$PRIVA_BINARY"
    # The native files are written below cache_root in the folders these keys
    # name; the steps' own folders lie outside it, so they are removed here.
    # vep_plugins_cachedir is restored below: the VEP step checks it exists.
    for k in clinvar_vcf_dir vep_plugins_cachedir gnomad_vcf_dir cadd_parent_dir loftee_parent_dir \
             hpo_raw_dir gene_mechanism_cache_dir gene_mechanism_raw_dir gofcards_workdir; do
        sed -i -E "/^${k}[[:space:]]*:/d" "$merged"
    done
    mkdir -p "$work/VEP_plugin_cache"
    if ((stub == 0)); then
        source "$SCRIPT_DIR/install_utils.sh" || die "cannot source install_utils.sh"
        export TMPDIR=$PWD/tmp
        "${PRIVA_BINARY:?PRIVA_BINARY must name the priva program}" cache-install prepare \
            --config "$merged" --build-config "$work/build.yaml" --source-root "$BASE_DIR" || die "cache-install prepare failed"
        "$PRIVA_BINARY" cache-install publish \
            --config "$merged" --build-config "$work/build.yaml" --source-root "$BASE_DIR" || die "cache-install publish failed"
    fi
    set_key "$merged" vep_plugins_cachedir "$work/VEP_plugin_cache"
    cp "$merged" config.yaml
    jq -n --arg l "$lineage" --arg done "$(date -Is)" --argjson steps "$steps" \
        --arg pi "$(jq -r .priva_identity "$plan")" \
        '{schema: 1, lineage: $l, published: $done, priva_identity: $pi, steps: $steps}' > resources.lock.json
}

cmd_adopt() {
    local plan id lineage receipt
    while (($#)); do case $1 in
        --plan) plan=$2; shift 2 ;; --step) id=$2; shift 2 ;; --lineage) lineage=$2; shift 2 ;;
        --receipt) receipt=$2; shift 2 ;; *) die "adopt: unknown option $1" ;; esac; done
    local a; a=$(jq -c --arg id "$id" '.steps[$id].adopt' "$plan")
    [[ $a != null ]] || die "$id is not adopted in $plan"
    jq -n --arg id "$id" --arg l "$lineage" --argjson a "$a" --arg done "$(date -Is)" \
        '{schema: 1, id: $id, lineage: $l, adopted: true, keys: $a.keys, files: $a.identity, inputs: [], finished: $done}' > "$receipt"
}

SCRIPT_DIR=$HERE
case ${1:-} in
    plan) shift; cmd_plan "$@" ;;
    run) shift; cmd_run "$@" ;;
    adopt) shift; cmd_adopt "$@" ;;
    publish) shift; cmd_publish "$@" ;;
    *) die "usage: install_task.sh plan|run|adopt|publish ..." ;;
esac
