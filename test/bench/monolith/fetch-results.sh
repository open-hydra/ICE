#!/usr/bin/env bash
# Bring a finished job back: its .out file and every per-run ICE log, into
# <dest>/<jobid>/. Usage: fetch-results.sh <jobid> [<dest>]
#
# The 2026-09 campaign kept the per-run logs only on the cluster, under an
# account nobody else can read; the banners and timing lines of its four jobs
# are lost to us. Every job of this campaign is fetched whole.
set -eu
JOB=${1:?job id}
REMOTE=${ICE_PERF_REMOTE:-monolith}
ROOT=${ICE_PERF_ROOT_REMOTE:-/data_fast/gpassarani/ice-perf}
DEST=${2:-/data10/passarani/Desktop/Software/ice-perf/runs}/$JOB
mkdir -p "$DEST/logs"

rsync -a "$REMOTE:$ROOT/logs/*-$JOB.out" "$DEST/"
out=$(ls "$DEST"/*-"$JOB".out | head -1)
name=$(basename "$out" | sed "s/-$JOB.out//")      # ice-final | ice-omp | ice-smoke
W=$ROOT/${name#ice-}
# Per-run logs are named <tag>-r<rep>.log with no job id: take the ones the
# .out file names, so two jobs sharing a work directory do not mix.
grep -o 'tag=[^ ]*' "$out" | sort -u | sed 's/tag=//' > "$DEST/tags.txt"
rsync -a --files-from=<(sed 's/$/-r*.log/' "$DEST/tags.txt" | sed "s#^#logs/#") \
  --no-implied-dirs "$REMOTE:$W/" "$DEST/" 2>/dev/null \
  || rsync -a "$REMOTE:$W/logs/" "$DEST/logs/"
echo "fetched $(ls "$DEST/logs" | wc -l) run logs and $(basename "$out") into $DEST"
