#!/usr/bin/env bash
#===============================================================================
#  Fast test: a probe given by position takes the cell whose centre is nearest.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
