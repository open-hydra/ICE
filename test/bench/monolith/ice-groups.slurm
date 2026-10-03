#!/bin/bash
#SBATCH --job-name=ice-groups
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Thread groups (ICE_THREAD_GROUPS): one rank of 80 threads with its blocks owned
# by groups of threads -- with OMP_PROC_BIND=close, one group per socket -- against
# the same binary with one group, the phase-2 binary (935af6b) and, as the target,
# the same cut run as 4 ranks x 20 (each rank's blocks on its own socket). Cuts of
# the 192^3 case: 24 blocks (24x64x192), 12, 4, and the single block (where the
# groups fall back to one: it must not move). Every run reports its windows, a
# steady time and a start-up percentage (_ice_common.sh).
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/groups
NEW=$ROOT/bin/ICE-mpi
REF=$ROOT/ref/ICE-mpi.935af6b
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
use_bin () {   # use_bin <binary>: snapshot it as this job's ICE-mpi and log its md5
  export ICE_PERF_BIN=$1
  source $ROOT/scripts/_ice_common.sh
  ice_setup
}
one () {   # one <tag> <binary> <groups> <set> <ranks> <threads>
  local tag=$1 bin=$2 g=$3 set=$4 R=$5 T=$6
  use_bin $bin
  export ICE_THREAD_GROUPS=$g
  if [ $R -eq 1 ]; then ice_run "$tag" $set 1 $T 1 "$(spread $T)" 2
  else ice_run "$tag" $set $R $T $R NONE 2; fi
  unset ICE_THREAD_GROUPS
}

echo "=== 24 blocks ==="
one r24-1x80-g4   $NEW 4 INPUT-R24 1 80
one r24-1x80-g1   $NEW 1 INPUT-R24 1 80
one r24-1x80-ref  $REF 1 INPUT-R24 1 80
one r24-4x20      $NEW 1 INPUT-R24 4 20
echo "=== 12 blocks ==="
one r12-1x80-g4   $NEW 4 INPUT-R12 1 80
one r12-1x80-g1   $NEW 1 INPUT-R12 1 80
echo "=== 4 blocks ==="
one r4-1x80-g4    $NEW 4 INPUT-R4  1 80
one r4-1x80-g1    $NEW 1 INPUT-R4  1 80
one r4-4x20       $NEW 1 INPUT-R4  4 20
echo "=== one block (groups fall back to one) ==="
one r1-1x80-g4    $NEW 4 INPUT     1 80
one r1-1x80-ref   $REF 1 INPUT     1 80
echo "ICE-GROUPS DONE"
