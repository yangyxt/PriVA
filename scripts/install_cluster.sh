#!/usr/bin/env bash
# Staged, noninteractive cluster deployment; does not run a family analysis.
set -euo pipefail

script_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
manager=${MAMBA_EXE:-/stornext/snfs170/next-gen/scratch/yangyxt/mamba/bin/micromamba}

usage() {
    cat <<'EOF'
Usage: install_cluster.sh environment STATE_DIRECTORY
       install_cluster.sh source SNAPSHOT_DIRECTORY [REFERENCE_COMMIT]

environment installs/checks the isolated priva_acmg environment and writes a
resolved explicit package lock and a validated software receipt. Existing
environments are checked, never silently upgraded. Run installation in Slurm.
source exports the committed reference source to an immutable-version directory;
its data/ cache builders may subsequently create assembly-specific artifacts.
EOF
}

case ${1:-} in
environment)
    [[ $# == 2 ]] || { usage >&2; exit 64; }
    state_dir=$2
    mkdir -p -- "$state_dir"
    exec 9>"$state_dir/environment.install.lock"
    flock -n 9 || { echo 'Another PriVA environment installation is active' >&2; exit 75; }
    [[ -x $manager ]] || { echo "Missing manager: $manager" >&2; exit 69; }
    env_spec="$script_root/envs/priva-cluster.yml"
    if ! "$manager" env list | awk '$1 == "priva_acmg" { found=1 } END { exit !found }'; then
        "$manager" env create --name priva_acmg --file "$env_spec" --yes --strict-channel-priority
    fi
    activation=$("$manager" shell activate --shell bash --name priva_acmg)
    eval "$activation"
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} MPLBACKEND=Agg
    "$manager" list --name priva_acmg --explicit > "$state_dir/environment.explicit.txt.partial"
    mv -- "$state_dir/environment.explicit.txt.partial" "$state_dir/environment.explicit.txt"
    python "$script_root/scripts/check_cluster_environment.py" --spec "$env_spec" --output "$state_dir/environment.ready.json"
    ;;
source)
    [[ $# == 2 || $# == 3 ]] || { usage >&2; exit 64; }
    destination=$2
    revision=${3:-09223ef5a8b0ceb4c6d61055c6eefe62b51ba903}
    revision=$(git -C "$script_root" rev-parse --verify "${revision}^{commit}")
    if [[ -e $destination ]]; then
        [[ -f $destination/REFERENCE_COMMIT && $(cat "$destination/REFERENCE_COMMIT") == "$revision" ]] || {
            echo "Refusing an unrelated existing source directory: $destination" >&2; exit 65;
        }
        echo "Reusing scientific reference $revision at $destination"
    else
        staging=$(mktemp -d "${destination}.partial.XXXXXXXX")
        git -C "$script_root" archive "$revision" | tar -xf - -C "$staging"
        printf '%s\n' "$revision" > "$staging/REFERENCE_COMMIT"
        mv -- "$staging" "$destination"
    fi
    ;;
-h|--help) usage ;;
*) usage >&2; exit 64 ;;
esac
