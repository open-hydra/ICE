#!/usr/bin/env bash
#===============================================================================
#  Fast test: components odd across mirror faces one cell apart stay zero.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
