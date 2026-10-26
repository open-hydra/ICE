#!/usr/bin/env bash
#===============================================================================
#  Fast test: res-threshold relative to the first non-zero density residual.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
