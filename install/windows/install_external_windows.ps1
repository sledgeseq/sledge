param(
  [Parameter(Mandatory = $true)]
  [string]$SledgeDir
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path -Path $SledgeDir -PathType Container)) {
  throw "SledgeDir does not exist: $SledgeDir"
}

$resolved = (Resolve-Path $SledgeDir).Path
$linuxInstaller = "$resolved/install/linux/install_external_linux.sh"

Write-Host "[install_external_windows] Windows support is provided via WSL2."
Write-Host "[install_external_windows] Repo path: $resolved"

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
  throw "wsl.exe not found. Install WSL2 and retry."
}

$wslRepo = (wsl.exe wslpath -a "$resolved").Trim()
$wslInstaller = "$wslRepo/install/linux/install_external_linux.sh"

Write-Host "[install_external_windows] Running Linux installer inside WSL:"
Write-Host "  bash $wslInstaller $wslRepo"

wsl.exe bash -lc "bash '$wslInstaller' '$wslRepo'"

Write-Host ""
Write-Host "[install_external_windows] Done."
Write-Host "Use the generated external_tools paths from inside WSL for sledge_filter config."
