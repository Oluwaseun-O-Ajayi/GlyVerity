#!/bin/bash
# Run thermo_channels.exe with the mono and DLLs of the ThermoRawFileParser environment.
# Usage: bash tools/thermo_channels/run.sh run.raw [output folder]
set -euo pipefail
ENV="${ENV:-$HOME/envs/trfp}"
HERE="$(cd "$(dirname "$0")" && pwd)"
MONO_PATH="$(cat "$HERE/dll_dir.txt")" "$ENV/bin/mono" "$HERE/thermo_channels.exe" "$@"
