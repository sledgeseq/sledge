param(
  [Parameter(Mandatory = $true)]
  [string]$SledgeDir
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path -Path $SledgeDir -PathType Container)) {
  throw "SledgeDir does not exist: $SledgeDir"
}

if (-not (Get-Command wsl.exe -ErrorAction SilentlyContinue)) {
  throw "wsl.exe not found. Install WSL2 and retry."
}

$resolved = (Resolve-Path $SledgeDir).Path
$wslRepo = (wsl.exe wslpath -a "$resolved").Trim()
$wslScript = "$wslRepo/install/test_installation.sh"

Write-Host "[test_installation_windows] Running install test in WSL2:"
Write-Host "  bash $wslScript $wslRepo"

wsl.exe bash -lc "bash '$wslScript' '$wslRepo'"

Write-Host ""
Write-Host "[test_installation_windows] Done."
