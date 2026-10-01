#!/usr/bin/env bash
#SBATCH --job-name=ice-sched2
#SBATCH --partition=free
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:30:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Second scheduling leg. ice-sched (303564) found that asking the flux kernel
# for four tiles per thread at 1 x 80 -- which on a 192-plane block switches it
# to the shallow tiling, plane pieces with every k-face precomputed once into
# the shared kface array and no private plane buffers -- cut the iteration from
# 0.374 to 0.304 s under a STATIC schedule, while dynamic,1 made the cell loops
# pay per cell. This job measures the shallow path under static scheduling
# where the deep path was good: one socket, four ranks of twenty, three nodes,
# and the 24-block cut, with 2, 4 and 8 tiles per thread; the one-tile static
# run of each is the reference inside the job.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/sched2
source $ROOT/scripts/_ice_common.sh
ice_setup

spread () { python3 - "$1" <<'PY'
import subprocess, collections, sys
N=int(sys.argv[1])
out=subprocess.run(['lscpu','-p=CPU,SOCKET'],stdout=subprocess.PIPE).stdout.decode()
sock=collections.OrderedDict()
for line in out.splitlines():
    if line.startswith('#') or not line.strip(): continue
    c,s=line.split(',')[:2]; sock.setdefault(int(s),[]).append(int(c))
for s in sock: sock[s].sort()
S=len(sock)
if N<S: cores=sorted(sock[list(sock)[i]][0] for i in range(N))
else:   cores=sorted(c for s in sock for c in sock[s][:N//S])
print(','.join(map(str,cores)))
PY
}
unset OMP_SCHEDULE      # static, as the loops were written

echo "=== 1 x 80: 1, 2, 4, 8 tiles per thread, static ==="
for tpt in 1 2 4 8; do
  export ICE_TILES_PER_THREAD=$tpt
  ice_run "s2-1x80-t$tpt" INPUT 1 80 1 "$(spread 80)" 2
done

echo "=== 1 x 20 on one socket: does the shallow path cost where the deep one was good? ==="
for tpt in 1 4 8; do
  export ICE_TILES_PER_THREAD=$tpt
  ice_run "s2-1x20-t$tpt" INPUT 1 20 1 "$(spread 20)" 2
done

echo "=== 4 x 20 (96 x 96 x 192 blocks): deep today ==="
for tpt in 1 4; do
  export ICE_TILES_PER_THREAD=$tpt
  ice_run "s2-4x20-t$tpt" INPUT-R4 4 20 4 NONE 2
done

echo "=== 12 x 20 on three nodes (48 x 64 x 192 blocks) ==="
for tpt in 1 4; do
  export ICE_TILES_PER_THREAD=$tpt
  ice_run "s2-12x20-t$tpt" INPUT-R12 12 20 4 NONE 2
done

echo "=== the 24-block cut on one rank ==="
for tpt in 1 4; do
  export ICE_TILES_PER_THREAD=$tpt
  ice_run "s2-1x80-R24-t$tpt" INPUT-R24 1 80 1 "$(spread 80)" 2
done

echo "ICE-SCHED2 DONE"
