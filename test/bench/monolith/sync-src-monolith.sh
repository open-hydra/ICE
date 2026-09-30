#!/usr/bin/env bash
# Push THIS ICE tree (the checkout this script lives in, submodules included) to
# monolith and rebuild the campaign binary there. Build artefacts stay on the
# cluster; only sources cross. Usage: sync-src-monolith.sh [cfg ...]  (default: mpi)
#
# The tree travels as files, not as a git push: the branch is pushed to origin
# only when the whole plan is finished and checked (user, 2026-09-30). The
# commit, the dirty-file count and the submodule pins are written to
# ICE-commit.txt so every job log can name what it ran.
#
# TecIO is a third-party C++ dependency built in place, per host, with that
# host's compilers -- the install directory is even named after them. Its
# SOURCES must cross (they live in ORION's tree); its build and install
# directories must not be deleted by --delete. Excluded paths are protected
# from --delete, so naming them is enough.
set -eu
REMOTE=${ICE_PERF_REMOTE:-monolith}
DEST=${ICE_PERF_ROOT_REMOTE:-/data_fast/gpassarani/ice-perf}
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=$(cd "$HERE/../../.." && pwd)          # the ICE root this script belongs to

{
  echo "ICE $(git -C "$SRC" rev-parse HEAD) dirty=$(git -C "$SRC" status --porcelain | wc -l)"
  git -C "$SRC" submodule status | awk '{print $2, $1}'
  echo "synced $(date -u +%Y-%m-%dT%H:%M:%SZ) from $(hostname):$SRC"
} > "$SRC/ICE-commit.txt"

ssh "$REMOTE" "mkdir -p $DEST/src $DEST/scripts $DEST/bin $DEST/logs"

rsync -a --delete \
  --exclude 'lib/ORION/lib/TecIO/tecio*-build-*' \
  --exclude 'lib/ORION/lib/TecIO/tecio*-install-*' \
  --exclude '.git' --exclude '.gitmodules' --exclude '.claude' \
  --exclude 'build*' --exclude 'bin' --exclude 'bin-*' --exclude 'modules' \
  --exclude '*.o' --exclude '*.mod' --exclude '*.a' \
  --exclude 'docs/' --exclude 'site/' --exclude 'plan-bucket/' \
  --exclude 'test/*/*/work' --exclude 'test/*/*/OUTPUT' \
  "$SRC/" "$REMOTE:$DEST/src/"

# The harness itself, so a job on the cluster runs the scripts of this commit.
rsync -a "$HERE/" "$REMOTE:$DEST/scripts/"

# hyb (MPI + OpenMP) is the campaign binary; the harness runs $DEST/bin/ICE-mpi.
for cfg in "${@:-hyb}"; do
  ssh "$REMOTE" "bash -l -c 'CC=gcc bash $DEST/scripts/build-cfg.sh $cfg $DEST/src'"
  ssh "$REMOTE" "cp $DEST/src/bin-$cfg/ICE $DEST/bin/ICE-$cfg && cp $DEST/src/ICE-commit.txt $DEST/bin/ICE-$cfg.commit \
                 && { [ $cfg = hyb ] && cp $DEST/bin/ICE-hyb $DEST/bin/ICE-mpi && cp $DEST/bin/ICE-hyb.commit $DEST/bin/ICE-mpi.commit; true; } \
                 && md5sum $DEST/bin/ICE-$cfg"
done
