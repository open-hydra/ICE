#!/usr/bin/env bash
#===============================================================================
#  Fast test: a restart goes on from the time of the solution it reads.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
