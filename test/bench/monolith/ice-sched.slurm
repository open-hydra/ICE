#!/usr/bin/env bash
#SBATCH --job-name=ice-sched
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:30:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The scheduling leg: the VTune profile of one rank of 80 threads (job 303562)
# put 24 % of its CPU time in spin -- the barriers of the static worksharing
# loops and the fork barrier between regions -- against 8 % with 20 threads.
# The step's loops are SCHEDULE(RUNTIME) since the instrument commit, so this
# job runs the same binary under several OMP_SCHEDULE values and tile counts
# (ICE_TILES_PER_THREAD) at 1 x 80, with 1 x 20 and 4 x 20 as references. Same
# protocol as the other jobs; the job's own anchor is not needed: the question
# is the ratio between legs of the same job.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/sched
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
packed () { python3 -c "print(','.join(str(i) for i in range($1)))"; }

leg () {   # leg <tag> <OMP_SCHEDULE or -> <tiles per thread> <threads> [set]
  local tag=$1 sched=$2 tpt=$3 T=$4 set=${5:-INPUT}
  if [ "$sched" = "-" ]; then unset OMP_SCHEDULE; else export OMP_SCHEDULE=$sched; fi
  export ICE_TILES_PER_THREAD=$tpt
  echo "SCHED tag=$tag OMP_SCHEDULE=${OMP_SCHEDULE:-unset} ICE_TILES_PER_THREAD=$tpt"
  ice_run "$tag" $set 1 $T 1 "$(spread $T)" 2
}

echo "=== 1 x 80: the schedule and the tile count ==="
leg sch-1x80-static-t1      -              1 80
leg sch-1x80-static-t4      static         4 80
leg sch-1x80-dyn1-t1        dynamic,1      1 80
leg sch-1x80-dyn1-t2        dynamic,1      2 80
leg sch-1x80-dyn1-t4        dynamic,1      4 80
leg sch-1x80-dyn4096-t4     dynamic,4096   4 80
leg sch-1x80-guided-t4      guided         4 80

echo "=== 1 x 20 on one socket: does a dynamic schedule cost where the static one was fine? ==="
leg sch-1x20-static-t1      -              1 20
leg sch-1x20-dyn1-t4        dynamic,1      4 20

echo "=== the 24-block cut on one rank ==="
leg sch-1x80-R24-static-t1  -              1 80 INPUT-R24
leg sch-1x80-R24-dyn1-t4    dynamic,1      4 80 INPUT-R24

# The per-cell phases of the 24-block cut run 2x the CPU time of the 4-block cut
# at 1 x 80 (VTune, job 303563) while the flux sweep does not. A plane of a
# 24x64x192 block is 86 kB: one 2 MB huge page spans ~23 planes, the k-ranges of
# ~10 threads on up to four sockets, so no placement makes it local; a 96x96
# plane is 516 kB. Twenty threads packed on ONE socket cannot pay that; spread
# over four they can.
echo "=== the 24-block cut on 20 threads: one socket against four ==="
cat /sys/kernel/mm/transparent_hugepage/enabled /sys/kernel/mm/transparent_hugepage/defrag 2>/dev/null | tr '\n' ' '; echo
unset OMP_SCHEDULE; export ICE_TILES_PER_THREAD=1
ice_run "sch-1x20-R24-packed"  INPUT-R24 1 20 1 "$(packed 20)" 2
ice_run "sch-1x20-R24-spread"  INPUT-R24 1 20 1 "$(spread 20)" 2
ice_run "sch-1x20-R4-packed"   INPUT     1 20 1 "$(packed 20)" 2

echo "=== 4 x 20 reference, static and dynamic ==="
unset OMP_SCHEDULE; export ICE_TILES_PER_THREAD=1
ice_run "sch-4x20-static-t1" INPUT-R4 4 20 4 NONE 2
export OMP_SCHEDULE=dynamic,1 ICE_TILES_PER_THREAD=4
ice_run "sch-4x20-dyn1-t4"   INPUT-R4 4 20 4 NONE 2

echo "ICE-SCHED DONE"
