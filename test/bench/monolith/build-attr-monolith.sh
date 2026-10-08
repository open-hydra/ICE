#!/usr/bin/env bash
# Build, on monolith, the hyb binaries of the commits ice-attr.slurm compares, from the
# trees sync-attr-monolith.sh shipped to $ROOT/attr/src-<commit> (no git here), all
# with this commit's build-cfg.sh hyb. Run on the login node through `bash -l`:
#   ssh monolith "bash -l -c 'bash /data_fast/gpassarani/ice-perf/scripts/build-attr-monolith.sh'"
# FiNeR v2.1.0 (the commits from 037e94f on) clones its five dependencies at configure
# time, unpinned, when they are missing; every such tree gets instead the set the
# campaign tree $ROOT/src was built with (the binaries of jobs 303832/303833), so no
# tree depends on GitHub's HEAD of the day. Binaries: $ROOT/attr/bin/ICE-mpi.<commit>.
set -eu
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
A=$ROOT/attr
COMMITS=${*:-935af6b 037e94f e69ebcc fddd1cb 055f6d4 f50efc2}
CFG=$ROOT/scripts/build-cfg.sh
DEPS=$ROOT/src/lib/third_party/FiNeR/src/third_party
mkdir -p $A/bin
for c in $COMMITS; do
  S=$A/src-$c
  [ -f $S/ICE-commit.txt ] || { echo "no tree $S (run sync-attr-monolith.sh)"; exit 1; }
  if grep -q 'ensure_finer_dependencies' $S/CMakeLists.txt; then
    for d in BeFoR64 FACE FLAP PENF StringiFor; do
      [ -f $S/lib/third_party/FiNeR/src/third_party/$d/CMakeLists.txt ] || \
        { mkdir -p $S/lib/third_party/FiNeR/src/third_party && cp -a $DEPS/$d $S/lib/third_party/FiNeR/src/third_party/; }
    done
  fi
  echo "== $c: $(head -1 $S/ICE-commit.txt)"
  CC=gcc bash $CFG hyb $S > $A/build-$c.log 2>&1 || { echo "BUILD FAILED $c (log $A/build-$c.log)"; tail -20 $A/build-$c.log; exit 1; }
  cp $S/bin-hyb/ICE $A/bin/ICE-mpi.$c
  { cat $S/ICE-commit.txt; echo "built $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname) with build-cfg.sh hyb"; } > $A/bin/ICE-mpi.$c.commit
  md5sum $A/bin/ICE-mpi.$c
done
