#!/usr/bin/env bash
#===============================================================================
#  Fast test: an inlet injects into the domain on a face whose normal points out of it.
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
