#!/bin/bash
#SBATCH --job-name=ice-reduced2
#SBATCH --partition=free
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The phase-2 verdict: the reduced matrix (ice-reduced.slurm's points) on the
# committed tree -- C3 (the parallel first touch, re-landed) and C5d (dynamic
# flux tiles over shallow pieces, cell loops static set by ICE) -- with the
# published "final" binary of job 303485 (tree 461eb4f = the wave plus C3)
# run beside it in the same job, on the same nodes, point by point. The
# harness now prints windows=, steady_iter and transient_pct for every run.
# Two stress legs that are not production layouts: the 24-block cut on one
# rank of 80 threads (C5d+C3 0.71 s against 0.63 for 461eb4f, job 303746),
# and the same under numactl --membind=0-3 (first-touch placement kept, the
# task taken out of automatic NUMA balancing): flat at ~0.44 means the
# balancer fights the plane-level placement of small blocks; ~0.71 means the
# placement itself is wrong for many small blocks.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/reduced2
NEW=$ROOT/bin/ICE-mpi
REF=$ROOT/final/bin.303485/ICE-mpi
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
packed () { python3 -c "print(','.join(str(i) for i in range($1)))"; }
use_bin () {   # use_bin <binary>: snapshot it as this job's ICE-mpi and log its md5
  export ICE_PERF_BIN=$1
  source $ROOT/scripts/_ice_common.sh
  ice_setup
}
pair () {   # pair <tag> <args of ice_run after the tag>: the new binary, then the reference
  local tag=$1; shift
  use_bin $NEW; unset ICE_WRAP; ice_run "new-$tag" "$@"
  use_bin $REF; unset ICE_WRAP; ice_run "ref-$tag" "$@"
}

echo "=== anchor ==="
pair env-C1-1x1          INPUT     1  1  1 "$(spread 1)"  2
echo "=== one node ==="
pair env-C80-4x20        INPUT-R4  4 20  4 "$(spread 80)" 3
pair omp-C80-1x80        INPUT     1 80  1 "$(spread 80)" 2
echo "=== three nodes ==="
pair env-C240-12x20      INPUT-R12 12 20 4 NONE 2
echo "=== one socket ==="
pair sock-C20-1x20-packed INPUT    1 20  1 "$(packed 20)" 2
echo "=== stress: 24 blocks on one rank of 80 threads ==="
pair blk-C80-1x80-R24    INPUT-R24 1 80  1 "$(spread 80)" 2
use_bin $NEW; export ICE_WRAP="numactl --membind=0-3"
ice_run "new-blk-C80-1x80-R24-membind" INPUT-R24 1 80 1 "$(spread 80)" 2
unset ICE_WRAP
echo "ICE-REDUCED2 DONE"
