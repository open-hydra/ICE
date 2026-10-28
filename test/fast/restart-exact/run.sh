#!/usr/bin/env bash
#===============================================================================
#  Fast test: a restart goes on as the uninterrupted run would.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
