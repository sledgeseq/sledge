#!/usr/bin/env bash
# Install MMseqs2, BLAST+, and FASTA36 for macOS.
# Uses Homebrew for MMseqs2/BLAST+ and builds FASTA36 from source.

set -euo pipefail

die() { echo "[install_external_macos] ERROR: $*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"; }

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <sledge_dir>" >&2
  exit 2
fi
[[ -d "$1" ]] || die "not a directory: $1"
SLEDGE_DIR="$(cd "$1" && pwd)"
[[ -f "${SLEDGE_DIR}/Makefile" ]] || die "does not look like sledge root: ${SLEDGE_DIR}"
[[ "$(uname -s)" == "Darwin" ]] || die "this script is macOS-only"

EXTERNAL="${SLEDGE_DIR}/external_tools"
WORKDIR="${EXTERNAL}/.downloads"
mkdir -p "${WORKDIR}"

need_cmd brew
need_cmd git
need_cmd make

FASTA36_REF="${FASTA36_REF:-master}"
FASTA36_REPO="https://github.com/wrpearson/fasta36.git"

mkdir -p "${EXTERNAL}/mmseqs/bin" "${EXTERNAL}/ncbi-blast/bin" "${EXTERNAL}/fasta36/bin"

if [[ "${SKIP_MMSEQS:-0}" != "1" ]]; then
  echo "[install_external_macos] Installing MMseqs2 via Homebrew..."
  brew list mmseqs2 >/dev/null 2>&1 || brew install mmseqs2
  mmseqs_path="$(command -v mmseqs || true)"
  [[ -x "${mmseqs_path}" ]] || die "mmseqs not found after brew install"
  ln -sf "${mmseqs_path}" "${EXTERNAL}/mmseqs/bin/mmseqs"
fi

if [[ "${SKIP_BLAST:-0}" != "1" ]]; then
  echo "[install_external_macos] Installing BLAST+ via Homebrew..."
  brew list blast >/dev/null 2>&1 || brew install blast
  blastp_path="$(command -v blastp || true)"
  makeblastdb_path="$(command -v makeblastdb || true)"
  [[ -x "${blastp_path}" && -x "${makeblastdb_path}" ]] || die "blast tools not found after brew install"
  ln -sf "${blastp_path}" "${EXTERNAL}/ncbi-blast/bin/blastp"
  ln -sf "${makeblastdb_path}" "${EXTERNAL}/ncbi-blast/bin/makeblastdb"
fi

if [[ "${SKIP_FASTA:-0}" != "1" ]]; then
  echo "[install_external_macos] Building FASTA36..."
  fasta_src="${WORKDIR}/fasta36"
  rm -rf "${fasta_src}" "${EXTERNAL}/fasta36-src"
  git clone --depth 1 --branch "${FASTA36_REF}" "${FASTA36_REPO}" "${fasta_src}" 2>/dev/null || git clone --depth 1 "${FASTA36_REPO}" "${fasta_src}"
  mv "${fasta_src}" "${EXTERNAL}/fasta36-src"

  pushd "${EXTERNAL}/fasta36-src/src" >/dev/null
  built=0
  for mk in ../make/Makefile.os_x86_64 ../make/Makefile.os_x86_64_sse4 ../make/Makefile.os_darwin ../make/Makefile.macosx; do
    if [[ -f "${mk}" ]]; then
      make -f "${mk}" -j"$(sysctl -n hw.ncpu 2>/dev/null || echo 4)" all && built=1 && break
    fi
  done
  popd >/dev/null
  [[ "${built}" -eq 1 ]] || die "fasta36 build failed for macOS"
  [[ -x "${EXTERNAL}/fasta36-src/bin/ssearch36" ]] || die "ssearch36 missing after build"
  ln -sf "${EXTERNAL}/fasta36-src/bin/ssearch36" "${EXTERNAL}/fasta36/bin/ssearch36"
fi

echo ""
echo "Done. Point sledge_filter config at:"
echo "  MMSEQS=\"${EXTERNAL}/mmseqs/bin/mmseqs\""
echo "  BLAST_DIR=\"${EXTERNAL}/ncbi-blast/bin\""
echo "  FASTA_DIR=\"${EXTERNAL}/fasta36/bin\""
