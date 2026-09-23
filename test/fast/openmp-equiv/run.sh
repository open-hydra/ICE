#!/usr/bin/env bash
#===============================================================================
#  Fast test: the number of OpenMP threads must not change the result.
#
#  Runs a short Doisneau/IG (crossing jets, IG closure) with 1 and with 4
#  threads and requires the two solutions to be bit-identical. A data race in
#  the flux or boundary loops shows up here as a difference.
#===============================================================================
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
source "$HERE/../common.sh"

WORK=$HERE/work
ice_fast_prepare "$WORK" "${ICE_FAST_CASE:-Doisneau/IG}" "${ICE_FAST_ITER:-200}" 1 4

for n in 1 4; do
  ice_fast_run "$WORK/run$n" "$n" 1
done

ice_fast_compare "$WORK/run1" "$WORK/run4" "1 thread" "4 threads"
