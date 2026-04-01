#!/usr/bin/env bash
# Dispatcher for OS-specific external tool installers.
#
# Usage:
#   install/install_external.sh <sledge_dir>
#
# Routes to:
#   install/linux/install_external_linux.sh
#   install/macos/install_external_macos.sh
#   install/windows/install_external_windows.ps1

set -euo pipefail

die() { echo "[install_external] ERROR: $*" >&2; exit 1; }

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <sledge_dir>" >&2
  exit 2
fi

[[ -d "$1" ]] || die "not a directory: $1"
SLEDGE_DIR="$(cd "$1" && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OS_NAME="$(uname -s)"

case "${OS_NAME}" in
  Linux)
    exec "${SCRIPT_DIR}/linux/install_external_linux.sh" "${SLEDGE_DIR}"
    ;;
  Darwin)
    exec "${SCRIPT_DIR}/macos/install_external_macos.sh" "${SLEDGE_DIR}"
    ;;
  MINGW*|MSYS*|CYGWIN*)
    if command -v powershell.exe >/dev/null 2>&1; then
      exec powershell.exe -ExecutionPolicy Bypass -File "${SCRIPT_DIR}/windows/install_external_windows.ps1" -SledgeDir "${SLEDGE_DIR}"
    fi
    die "Windows shell detected. Run install/windows/install_external_windows.ps1 from PowerShell."
    ;;
  *)
    die "unsupported OS: ${OS_NAME}"
    ;;
esac
