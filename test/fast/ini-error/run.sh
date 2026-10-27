#!/usr/bin/env bash
#===============================================================================
#  Fast test: input.ini edited during a run with an error, ini-error = warn | stop.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
