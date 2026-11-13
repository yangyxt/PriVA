#!/usr/bin/env bash
# Install/check the pinned driver only; no workflow or shell-startup changes.
# Downloads stay in a pin-bound sibling cache. Interrupted builds are retained;
# a retry builds afresh from verified bytes and publishes by directory rename.
set -euo pipefail

nf_url='https://github.com/nextflow-io/nextflow/releases/download/v25.10.4/nextflow-25.10.4-dist'
nf_sha='60db869a7618783e3ee60183a2b0c6c50ffc67221eabec60800b904991786aaa'
java_url='https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.12.1%2B1/OpenJDK21U-jre_x64_linux_hotspot_21.0.12.1_1.tar.gz'
java_sha='2413149700df0f7d440500a84a8f764c535f21e5a5e87d38328b64eec2c5b500'
java_directory='jdk-21.0.12.1+1-jre'
nf_filename='nextflow-25.10.4-dist'

die() { printf 'Nextflow runtime: %s\n' "$*" >&2; exit 1; }
usage() {
    printf 'Usage: bash %s install|check --prefix /absolute/runtime-prefix\n' "$0"
    printf 'Legacy compatibility: bash %s /absolute/runtime-prefix\n' "$0"
}

parse_prefix() {
    local parsed component current= option
    parsed=$(getopt -o '' --long prefix:,help -- "$@") || exit 64
    eval set -- "$parsed"
    prefix=
    while [[ $# -gt 0 ]]; do
        option=$1; shift
        case $option in
            --prefix) [[ -z $prefix ]] || die 'Supply --prefix once'; prefix=$1; shift ;;
            --help) usage; exit 0 ;;
            --) break ;;
        esac
    done
    [[ $# == 0 && $prefix == /* && $prefix != / && $prefix != *$'\n'* ]] || { usage >&2; exit 64; }
    # Check the supplied path before normalization, including existing ancestors.
    local -a components
    IFS=/ read -r -a components <<< "$prefix"
    for component in "${components[@]}"; do
        [[ -z $component || $component == . ]] && continue
        [[ $component != .. ]] || die 'Prefix must not contain .. components'
        current+=/$component
        [[ ! -L $current ]] || die "Refusing symlinked prefix component: $current"
    done
    prefix=$(realpath -m -- "$prefix")
    [[ $prefix != / ]] || die 'Refusing the filesystem root as a runtime prefix'
    [[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || die 'This lock supports Linux x86_64 only'
}

require_tools() {
    local program
    for program in "$@"; do command -v "$program" >/dev/null 2>&1 || die "Missing prerequisite: $program"; done
}
regular_file() { [[ -f $1 && ! -L $1 ]] || die "Expected a regular file, retained: $1"; }
physical_directory() { [[ -d $1 && ! -L $1 ]] || die "Expected a physical directory, retained: $1"; }
matches_hash() { [[ $(sha256sum -- "$1" | cut -d ' ' -f 1) == "$2" ]]; }
owner_marker() { printf '%s\n' 'PriVA Nextflow runtime v1' "$prefix" "$nf_url" "$nf_sha" "$java_url" "$java_sha"; }
checksums() { printf '%s  %s\n' "$nf_sha" "bin/$nf_filename" "$java_sha" 'downloads/temurin-jre.tar.gz'; }

launcher() {
    cat <<'LAUNCHER'
#!/usr/bin/env bash
set -euo pipefail
runtime=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
export JAVA_HOME="$runtime/jdk-21.0.12.1+1-jre"
export NXF_HOME="${NXF_HOME:-$runtime/nxf-home}"
export NXF_VER=25.10.4
export NXF_SYNTAX_PARSER=v2
export NXF_OPTS="${NXF_OPTS:--Xms128m -Xmx1536m -XX:ActiveProcessorCount=2}"
exec "$runtime/bin/nextflow-25.10.4-dist" "$@"
LAUNCHER
}

verify_jre() {
    # Compare extracted bytes/links with the checksum-verified tar: version text
    # alone cannot identify an altered JVM/library. Python3 is supplied by the
    # Conda prerequisite used by main_install; no third-party modules are used.
    python3 - "$1" "$java_directory" <<'PY'
import hashlib
import os
from pathlib import Path, PurePosixPath
import sys
import tarfile

root = Path(sys.argv[1])
jre = root / sys.argv[2]
expected = set()
def require(condition, message):
    if not condition:
        raise RuntimeError(message)
def digest(stream):
    value = hashlib.sha256()
    for block in iter(lambda: stream.read(1048576), b''):
        value.update(block)
    return value.digest()
try:
    with tarfile.open(root / 'downloads/temurin-jre.tar.gz', 'r:gz') as archive:
        for member in archive:
            relative = PurePosixPath(member.name)
            require(not relative.is_absolute() and '..' not in relative.parts
                    and relative.parts[0] == jre.name, 'Unexpected JRE archive path')
            path = root.joinpath(*relative.parts)
            expected.add(path)
            for parent in path.parents:
                if parent == root:
                    break
                require(not parent.is_symlink(), 'Symlinked JRE parent: ' + str(parent))
            if member.isdir():
                require(path.is_dir() and not path.is_symlink(), 'Missing JRE directory: ' + str(path))
            elif member.issym():
                target = Path(os.path.abspath(path.parent / member.linkname))
                require(os.path.commonpath((str(target), str(jre))) == str(jre)
                        and path.is_symlink() and os.readlink(path) == member.linkname,
                        'JRE symlink differs or escapes its root: ' + str(path))
            else:
                require(member.isfile() or member.islnk(), 'Unsupported JRE archive member')
                require(path.is_file() and not path.is_symlink(), 'Missing JRE file: ' + str(path))
                with archive.extractfile(member) as original, path.open('rb') as installed:
                    require(digest(original) == digest(installed), 'JRE bytes differ: ' + str(path))
    actual = {jre}
    for directory, dirs, files in os.walk(jre, followlinks=False):
        actual.update(Path(directory) / name for name in dirs + files)
    require(expected == actual, 'Unexpected or missing extracted JRE members')
except (OSError, ValueError, RuntimeError, tarfile.TarError) as error:
    sys.exit('Nextflow runtime: ' + str(error))
PY
}

verify_runtime() {
    local root=$1 output child probe_tmp
    physical_directory "$root"
    for child in bin downloads nxf-home "$java_directory"; do physical_directory "$root/$child"; done
    for child in bin/nextflow "bin/$nf_filename" downloads/temurin-jre.tar.gz SHA256SUMS; do regular_file "$root/$child"; done
    if [[ -e $root/.priva-nextflow-runtime || -L $root/.priva-nextflow-runtime ]]; then
        regular_file "$root/.priva-nextflow-runtime"
        cmp -s "$root/.priva-nextflow-runtime" <(owner_marker) || die 'Runtime ownership/pins differ; retained'
    fi
    cmp -s "$root/SHA256SUMS" <(checksums) || die 'Runtime checksum manifest differs; retained'
    matches_hash "$root/bin/$nf_filename" "$nf_sha" || die 'Nextflow distribution checksum differs; retained'
    matches_hash "$root/downloads/temurin-jre.tar.gz" "$java_sha" || die 'JRE archive checksum differs; retained'
    cmp -s "$root/bin/nextflow" <(launcher) || die 'Launcher differs; retained'
    [[ -x $root/bin/nextflow && -x $root/bin/$nf_filename && -x $root/$java_directory/bin/java ]] || die 'Runtime executable permission is missing'
    verify_jre "$root"
    [[ -z $(find "$root/nxf-home" -type l -print -quit) ]] || die 'Refusing symlinked Nextflow state'
    probe_tmp=${TMPDIR:-$root/nxf-home}
    output=$(timeout --kill-after=10s 90s env -u JAVA_TOOL_OPTIONS -u _JAVA_OPTIONS -u JDK_JAVA_OPTIONS \
        TMPDIR="$probe_tmp" "$root/$java_directory/bin/java" -Xms32m -Xmx256m -XX:ActiveProcessorCount=1 \
        "-Djava.io.tmpdir=$probe_tmp" -version 2>&1) || die 'Java version check failed'
    grep -Fq 'version "21.0.12.1"' <<< "$output" && grep -Fq 'Temurin-21.0.12.1+1' <<< "$output" || die 'Java version differs'
    # Persistent state stays with the runtime; caller-selected scratch must stay
    # in TMPDIR. The packaged launcher maps NXF_TEMP to JVM java.io.tmpdir.
    output=$(cd -- "$root/nxf-home" && timeout --kill-after=10s 90s env -u JAVA_TOOL_OPTIONS -u _JAVA_OPTIONS -u JDK_JAVA_OPTIONS \
        JAVA_HOME="$root/$java_directory" NXF_JAVA_HOME="$root/$java_directory" \
        NXF_HOME="$root/nxf-home" NXF_OFFLINE=true TMPDIR="$probe_tmp" NXF_TEMP="$probe_tmp" \
        XDG_CACHE_HOME="$root/nxf-home" NXF_OPTS='-Xms64m -Xmx768m -XX:ActiveProcessorCount=1' \
        "$root/bin/nextflow" -version 2>&1) || die 'Nextflow version check failed'
    grep -Eq 'version[[:space:]]+25[.]10[.]4([^0-9.]|$)' <<< "$output" || die 'Nextflow version differs'
}

prepare_download_cache() {
    local seed child
    if [[ ! -e $stage && ! -L $stage ]]; then
        # Publish ownership with the directory. Interrupted seeds are retained
        # but never reused or traversed by a later attempt.
        seed=$(mktemp -d "${prefix}.seed.XXXXXXXX")
        owner_marker > "$seed/.priva-nextflow-runtime"
        mv -T --no-clobber -- "$seed" "$stage"
        [[ ! -e $seed ]] || die "Download cache appeared concurrently; retained $seed"
    fi
    physical_directory "$stage"
    regular_file "$stage/.priva-nextflow-runtime"
    [[ -O $stage ]] && cmp -s "$stage/.priva-nextflow-runtime" <(owner_marker) || die 'Unrelated download cache/ownership/pins; retained'
    shopt -s nullglob dotglob
    for child in "$stage"/*; do
        case ${child##*/} in
            .priva-nextflow-runtime|"$nf_filename"|"$nf_filename.complete"|temurin-jre.tar.gz|temurin-jre.tar.gz.complete) regular_file "$child" ;;
            *) die "Unrelated staging child, retained: $child" ;;
        esac
    done
    shopt -u nullglob dotglob
}

download() {
    local destination=$1 url=$2 expected=$3
    if [[ -e $destination || -L $destination ]]; then
        regular_file "$destination"
        matches_hash "$destination" "$expected" && return 0
        [[ ! -e $destination.complete ]] || die "Completed download checksum differs; retained: $destination"
    fi
    # Resume interrupted owned bytes. A completed bad download is retained and
    # fails instead of being silently deleted or replaced.
    curl --fail --location --retry 3 --connect-timeout 30 --speed-limit 1024 --speed-time 120 \
        --continue-at - "$url" -o "$destination"
    printf '%s\n' "$expected" > "$destination.complete"
    matches_hash "$destination" "$expected" || die "Downloaded checksum differs; retained: $destination"
}

install() {
    parse_prefix "$@"
    require_tools curl sha256sum cut tar flock python3 cmp find grep env timeout mktemp ln mv
    local stage="${prefix}.install" build lock="${prefix}.install.lock"
    mkdir -p -- "$(dirname -- "$prefix")"
    [[ ! -L $lock && ( ! -e $lock || -f $lock ) ]] || die 'Refusing unrelated installation lock'
    exec 9>>"$lock"
    flock -n 9 || { printf 'Nextflow runtime: another installer holds %s\n' "$lock" >&2; exit 75; }
    if [[ -e $prefix || -L $prefix ]]; then
        verify_runtime "$prefix"
        printf 'Verified existing Nextflow runtime: %s\n' "$prefix"
        return
    fi
    prepare_download_cache
    download "$stage/$nf_filename" "$nf_url" "$nf_sha"
    download "$stage/temurin-jre.tar.gz" "$java_url" "$java_sha"
    build=$(mktemp -d "${prefix}.partial.XXXXXXXX")
    owner_marker > "$build/.priva-nextflow-runtime"
    mkdir -- "$build/bin" "$build/downloads" "$build/nxf-home"
    ln -- "$stage/$nf_filename" "$build/bin/$nf_filename"
    ln -- "$stage/temurin-jre.tar.gz" "$build/downloads/temurin-jre.tar.gz"
    tar -xzf "$build/downloads/temurin-jre.tar.gz" --no-same-owner -C "$build"
    launcher > "$build/bin/nextflow"
    chmod +x "$build/bin/nextflow" "$build/bin/$nf_filename"
    checksums > "$build/SHA256SUMS"
    verify_runtime "$build"
    mv -T --no-clobber -- "$build" "$prefix"
    [[ ! -e $build ]] || die "Prefix appeared before publication; complete build retained: $build"
    printf 'Installed pinned Nextflow runtime: %s\n' "$prefix"
}

check() {
    parse_prefix "$@"
    require_tools sha256sum cut python3 cmp find grep env timeout
    verify_runtime "$prefix"
    printf 'Verified pinned Nextflow runtime: %s\n' "$prefix"
}

main() {
    if [[ $# == 1 && $1 == /* ]]; then install --prefix "$1"; return; fi
    local procedure=${1:-}; [[ $# == 0 ]] || shift
    case $procedure in
        install|check) "$procedure" "$@" ;;
        -h|--help) usage ;;
        *) usage >&2; exit 64 ;;
    esac
}
main "$@"
