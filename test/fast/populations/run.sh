#!/usr/bin/env bash
#===============================================================================
#  Fast test: which copy of the boundary table belongs to which particle family.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
