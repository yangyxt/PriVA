#!/usr/bin/env bash
# Install the pinned Linux x86_64 workflow driver; no shell startup changes.
set -euo pipefail
if [[ $# != 1 || $1 != /* ]]; then
    echo "Usage: $0 /absolute/new/runtime-prefix" >&2
    exit 64
fi
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]] || {
    echo 'This lock supports Linux x86_64 only.' >&2; exit 64;
}
prefix=$1
[[ ! -e "$prefix" ]] || { echo "Refusing existing prefix: $prefix" >&2; exit 65; }
mkdir -p "$prefix/downloads" "$prefix/bin" "$prefix/nxf-home"
nf_url='https://github.com/nextflow-io/nextflow/releases/download/v25.10.4/nextflow-25.10.4-dist'
nf_sha='60db869a7618783e3ee60183a2b0c6c50ffc67221eabec60800b904991786aaa'
java_url='https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.12.1%2B1/OpenJDK21U-jre_x64_linux_hotspot_21.0.12.1_1.tar.gz'
java_sha='2413149700df0f7d440500a84a8f764c535f21e5a5e87d38328b64eec2c5b500'
curl --fail --location --retry 3 "$nf_url" -o "$prefix/bin/nextflow-25.10.4-dist"
curl --fail --location --retry 3 "$java_url" -o "$prefix/downloads/temurin-jre.tar.gz"
printf '%s  %s\n' "$nf_sha" "$prefix/bin/nextflow-25.10.4-dist" "$java_sha" "$prefix/downloads/temurin-jre.tar.gz" | sha256sum --check
tar -xzf "$prefix/downloads/temurin-jre.tar.gz" -C "$prefix"
chmod +x "$prefix/bin/nextflow-25.10.4-dist"
cat > "$prefix/bin/nextflow" <<'LAUNCHER'
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
chmod +x "$prefix/bin/nextflow"
printf '%s  %s\n' "$nf_sha" 'bin/nextflow-25.10.4-dist' "$java_sha" 'downloads/temurin-jre.tar.gz' > "$prefix/SHA256SUMS"
"$prefix/bin/nextflow" -version
