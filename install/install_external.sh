#!/bin/bash
# Download / build MMseqs2, NCBI BLAST+, and FASTA36 (ssearch36) into <sledge_dir>/external_tools/
#
# Usage: install_external.sh <sledge_dir>
#   sledge_dir  Root of the sledge tree (relative or absolute; stored internally as absolute).
#
# Optional environment overrides:
#   MMSEQS_TAG=18-8cc5c          MMseqs2 release tag (GitHub)
#   MMSEQS_ARCH=sse2|sse41|avx2  Prebuilt binary variant (default: sse2, widest CPU support)
#   BLAST_VERSION=2.15.0         NCBI BLAST+ version directory on FTP (no "+" in path)
#   FASTA36_REF=master           Git ref for wrpearson/fasta36
#   SKIP_MMSEQS=1, SKIP_BLAST=1, SKIP_FASTA=1  Skip one or more steps

set -euo pipefail

die() { echo "[install_external] ERROR: $*" >&2; exit 1; }

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <sledge_dir>" >&2
  exit 2
fi
[[ -d "$1" ]] || die "not a directory: $1"
SLEDGE_DIR="$(cd "$1" && pwd)"
[[ -f "${SLEDGE_DIR}/Makefile" || -f "${SLEDGE_DIR}/install/install_external.sh" ]] \
  || die "does not look like sledge root: ${SLEDGE_DIR}"

# Must be absolute: after `cd` into fasta36/src, relative paths like ../external_tools break checks below.
EXTERNAL="${SLEDGE_DIR}/external_tools"
WORKDIR="${EXTERNAL}/.downloads"
mkdir -p "${WORKDIR}"

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

fetch() {
  local url="$1" dest="$2"
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL --retry 3 --connect-timeout 30 -o "$dest" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget -q --tries=3 --timeout=30 -O "$dest" "$url"
  else
    die "need curl or wget to download"
  fi
}

uname_m="$(uname -m)"
[[ "${uname_m}" == "x86_64" ]] || die "this script targets Linux x86_64; detected: ${uname_m}"

need_cmd tar
need_cmd git
need_cmd make

MMSEQS_TAG="${MMSEQS_TAG:-18-8cc5c}"
MMSEQS_ARCH="${MMSEQS_ARCH:-sse2}"
BLAST_VERSION="${BLAST_VERSION:-2.15.0}"
FASTA36_REF="${FASTA36_REF:-master}"

MMSEQS_URL="https://github.com/soedinglab/MMseqs2/releases/download/${MMSEQS_TAG}/mmseqs-linux-${MMSEQS_ARCH}.tar.gz"
BLAST_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
FASTA36_REPO="https://github.com/wrpearson/fasta36.git"

echo "[install_external] Sledge root: ${SLEDGE_DIR}"
echo "[install_external] Install prefix: ${EXTERNAL}"
mkdir -p "${EXTERNAL}"

# --- MMseqs2 ---
if [[ "${SKIP_MMSEQS:-0}" != "1" ]]; then
  echo "[install_external] Installing MMseqs2 (${MMSEQS_TAG}, ${MMSEQS_ARCH})..."
  mmseqs_tgz="${WORKDIR}/mmseqs-linux-${MMSEQS_ARCH}.tar.gz"
  fetch "${MMSEQS_URL}" "${mmseqs_tgz}"
  rm -rf "${EXTERNAL}/mmseqs"
  tar xzf "${mmseqs_tgz}" -C "${EXTERNAL}"
  [[ -x "${EXTERNAL}/mmseqs/bin/mmseqs" ]] || die "MMseqs2 binary missing after extract"
  echo "[install_external] MMSEQS=${EXTERNAL}/mmseqs/bin/mmseqs"
else
  echo "[install_external] Skipping MMseqs2 (SKIP_MMSEQS=1)"
fi

# --- BLAST+ ---
if [[ "${SKIP_BLAST:-0}" != "1" ]]; then
  echo "[install_external] Installing NCBI BLAST+ (${BLAST_VERSION})..."
  blast_tgz="${WORKDIR}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
  fetch "${BLAST_URL}" "${blast_tgz}"
  rm -rf "${EXTERNAL}/ncbi-blast"
  tar xzf "${blast_tgz}" -C "${EXTERNAL}"
  blast_dir="$(find "${EXTERNAL}" -maxdepth 1 -type d -name 'ncbi-blast-*+' | head -1)"
  [[ -n "${blast_dir}" ]] || die "could not find extracted ncbi-blast-*+ directory"
  mv "${blast_dir}" "${EXTERNAL}/ncbi-blast"
  [[ -x "${EXTERNAL}/ncbi-blast/bin/blastp" && -x "${EXTERNAL}/ncbi-blast/bin/makeblastdb" ]] \
    || die "BLAST+ binaries missing under ${EXTERNAL}/ncbi-blast/bin"
  echo "[install_external] BLAST_DIR=${EXTERNAL}/ncbi-blast/bin"
else
  echo "[install_external] Skipping BLAST+ (SKIP_BLAST=1)"
fi

# --- FASTA36 (ssearch36) ---
if [[ "${SKIP_FASTA:-0}" != "1" ]]; then
  echo "[install_external] Building FASTA36 (${FASTA36_REF})..."
  fasta_src="${WORKDIR}/fasta36"
  rm -rf "${fasta_src}"
  git clone --depth 1 --branch "${FASTA36_REF}" "${FASTA36_REPO}" "${fasta_src}" 2>/dev/null \
    || git clone --depth 1 "${FASTA36_REPO}" "${fasta_src}"
  rm -rf "${EXTERNAL}/fasta36"
  mv "${fasta_src}" "${EXTERNAL}/fasta36"
  cd "${EXTERNAL}/fasta36/src"
  _built=0
  for mk in ../make/Makefile.linux64_sse2 ../make/Makefile.linux64 ../make/Makefile.linux_sse2; do
    if [[ -f "${mk}" ]]; then
      echo "[install_external] make -f ${mk} ..."
      make -f "${mk}" -j"$(nproc 2>/dev/null || echo 4)" all && _built=1 && break
    fi
  done
  [[ "${_built}" -eq 1 ]] || die "fasta36 build failed (tried linux64_sse2 / linux64 / linux_sse2 makefiles)"
  [[ -x "${EXTERNAL}/fasta36/bin/ssearch36" ]] \
    || die "ssearch36 missing or not executable after build: ${EXTERNAL}/fasta36/bin/ssearch36"
  echo "[install_external] FASTA_DIR=${EXTERNAL}/fasta36/bin"
else
  echo "[install_external] Skipping FASTA36 (SKIP_FASTA=1)"
fi

echo ""
echo "Done. Point sledge_filter configs at:"
echo "  MMSEQS=\"${EXTERNAL}/mmseqs/bin/mmseqs\""
echo "  BLAST_DIR=\"${EXTERNAL}/ncbi-blast/bin\""
echo "  FASTA_DIR=\"${EXTERNAL}/fasta36/bin\""
