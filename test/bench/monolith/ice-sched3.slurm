#!/usr/bin/env bash
#SBATCH --job-name=ice-sched3
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Third scheduling leg: the shallow tiling where the deep one is good. On the
# 192-plane block the shallow path is taken only when the tiles asked for
# exceed the plane count, i.e. from 10 tiles per thread at 20 threads (ice-sched2
# measured 4 and 8 there, both still deep: finer deep tiles, +14 %). One socket
# of 20 and four ranks of 20, deep (1 tile) against shallow (10 tiles), static.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/sched3
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
unset OMP_SCHEDULE
for tpt in 1 10; do
  export ICE_TILES_PER_THREAD=$tpt
  ice_run "s3-1x20-t$tpt" INPUT    1 20 1 "$(spread 20)" 2
  ice_run "s3-4x20-t$tpt" INPUT-R4 4 20 4 NONE 2
done
echo "ICE-SCHED3 DONE"
