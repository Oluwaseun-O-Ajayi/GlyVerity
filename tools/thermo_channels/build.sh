#!/bin/bash
# Build thermo_channels.exe against the RawFileReader DLLs of a ThermoRawFileParser conda
# environment. Usage: bash tools/thermo_channels/build.sh [env prefix, default ~/envs/trfp]
set -euo pipefail
ENV="${1:-$HOME/envs/trfp}"
DLL_DIR="$(dirname "$(find "$ENV" -name ThermoFisher.CommonCore.RawFileReader.dll | head -1)")"
if [ -z "$DLL_DIR" ] || [ "$DLL_DIR" = "." ]; then
  echo "RawFileReader DLLs not found under $ENV" >&2; exit 1
fi
HERE="$(cd "$(dirname "$0")" && pwd)"
"$ENV/bin/mcs" -r:"$DLL_DIR/ThermoFisher.CommonCore.Data.dll" \
  -r:"$DLL_DIR/ThermoFisher.CommonCore.RawFileReader.dll" \
  "$HERE/thermo_channels.cs" -out:"$HERE/thermo_channels.exe"
echo "$DLL_DIR" > "$HERE/dll_dir.txt"
echo "built $HERE/thermo_channels.exe (DLLs in $DLL_DIR)"
