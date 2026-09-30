#!/usr/bin/env bash
# Push the ICE scaling case to monolith (once; it is ~5 GB of szplt).
#
# The case is Marco Grossi's 192^3 set with its MDB cuts, copied from
# /data10/grossi/tmp/scaling-ice/case/build on 2026-09-30 (README.md). /scratch
# on monolith is full, so it goes to /data_fast. Incremental and checksummed:
# a second run moves nothing and proves the copy.
set -eu
REMOTE=${ICE_PERF_REMOTE:-monolith}
LOCAL=${ICE_PERF_CASE_LOCAL:-/data10/passarani/Desktop/Software/ice-perf/case/build}
DEST=${ICE_PERF_ROOT_REMOTE:-/data_fast/gpassarani/ice-perf}

ssh "$REMOTE" "mkdir -p $DEST/case"

echo "== case (this is the slow part) =="
rsync -a --checksum --info=progress2 --exclude '.try*' --exclude 'mdb*' \
  "$LOCAL"/INPUT "$LOCAL"/INPUT-R* \
  "$REMOTE:$DEST/case/"

echo "== check =="
ssh "$REMOTE" "cd $DEST && du -sh case && ls case && df -h /data_fast | tail -1"
