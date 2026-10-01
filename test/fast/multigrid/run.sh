#!/usr/bin/env bash
#===============================================================================
#  Fast test: grid sequencing runs, and neither threads nor ranks change it.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
