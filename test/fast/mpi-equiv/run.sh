#!/usr/bin/env bash
#===============================================================================
#  Fast test: the number of MPI ranks must not change the result.
#
#  Runs a short multi-block case on 1 and on 2 ranks and requires the two
#  solutions to be bit-identical. The case comes from ICE_FAST_CASE, so the
#  same script covers the connection (IG-split4) and the chimera (IG-chimera)
#  halo paths. Needs an MPI build.
#===============================================================================
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
source "$HERE/../common.sh"

WORK=$HERE/work
ice_fast_prepare "$WORK" "${ICE_FAST_CASE:-Doisneau/IG-split4}" "${ICE_FAST_ITER:-200}" 1 2

for n in 1 2; do
  ice_fast_run "$WORK/run$n" 2 "$n"
done

# Both builds write bin/ICE, so make sure the binary really is the MPI one:
# without MPI, mpirun just starts two independent serial solvers.
if ! grep -q "Number of ranks   -->    2" "$WORK/run2/log"; then
  echo "[fast] FAIL: $ICE_BIN did not run as 2 MPI ranks -- is it the MPI build?"
  exit 1
fi

ice_fast_compare "$WORK/run1" "$WORK/run2" "1 rank" "2 ranks"
