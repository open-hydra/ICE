#!/usr/bin/env bash
# Ship the source trees of the commits ice-attr.slurm compares to monolith, then build
# them there (build-attr-monolith.sh). Usage: sync-attr-monolith.sh [commit ...]
#
# The trees are made HERE, one git worktree per commit with its submodules, because
# monolith's https access to GitHub drops large clones (2026-10-08: FiNeR's nested
# PENF clone failed with "curl 18 ... early EOF" after 16 min). Only files cross,
# as in sync-src-monolith.sh; each tree carries ICE-commit.txt naming what it is.
# The worktrees stay registered in this repository: afterwards remove each with
# `git worktree remove --force $ICE_ATTR_TMP/src-<commit>`.
set -eu
REMOTE=${ICE_PERF_REMOTE:-monolith}
DEST=${ICE_PERF_ROOT_REMOTE:-/data_fast/gpassarani/ice-perf}
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=$(cd "$HERE/../../.." && pwd)
TMP=${ICE_ATTR_TMP:?ICE_ATTR_TMP: a scratch directory for the worktrees}
COMMITS=${*:-935af6b 037e94f e69ebcc fddd1cb 055f6d4 f50efc2}
mkdir -p "$TMP"
ssh "$REMOTE" "mkdir -p $DEST/attr/bin $DEST/scripts"
rsync -a "$HERE/" "$REMOTE:$DEST/scripts/"
for c in $COMMITS; do
  T=$TMP/src-$c
  if [ ! -d "$T" ]; then
    git -C "$SRC" worktree add -q --detach "$T" "$c"
    git -C "$T" submodule update -q --init --recursive
  fi
  { echo "ICE $(git -C "$T" rev-parse HEAD) dirty=$(git -C "$T" status --porcelain | wc -l)"
    git -C "$T" submodule status | awk '{print $2, $1}'
    echo "synced $(date -u +%Y-%m-%dT%H:%M:%SZ) from $(hostname):$T"; } > "$T/ICE-commit.txt"
  rsync -a --delete --exclude '.git' --exclude 'build*' --exclude 'bin' --exclude 'bin-*' \
    --exclude 'docs/' --exclude 'plan-bucket/' --exclude 'test/*/*/work' --exclude 'test/*/*/OUTPUT' \
    "$T/" "$REMOTE:$DEST/attr/src-$c/"
  echo "synced $c"
done
