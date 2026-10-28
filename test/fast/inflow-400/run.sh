#!/usr/bin/env bash
#===============================================================================
#  Fast test: the run says where particles enter through extrapolation faces (400).
#===============================================================================
set -euo pipefail
cd "$(dirname "$0")"
exec python3 -B run.py
