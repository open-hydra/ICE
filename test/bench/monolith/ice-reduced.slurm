#!/bin/bash
#SBATCH --job-name=ice-reduced
#SBATCH --partition=free
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The per-phase matrix: the points that decide whether a change bought anything,
# at a fifth of the node-hours of ice-final.slurm. One anchor, the full node in
# the two layouts that bracket the rank/thread split, the three-node point that
# shows the halo, and the 20-core socket pair. Same protocol as the full matrix:
# one job, one node set, the job's own 1-core anchor, minimum over repetitions,
# core-cycles beside seconds. Run the full matrix at the end of each wave.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/reduced
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

echo "=== anchor ==="
ice_run "env-C1-1x1"          INPUT     1  1  1 "$(spread 1)"  2

echo "=== the full node, two layouts ==="
ice_run "env-C80-4x20"        INPUT-R4  4 20  4 "$(spread 80)" 3
ice_run "omp-C80-1x80"        INPUT     1 80  1 "$(spread 80)" 2

echo "=== three nodes: the halo ==="
ice_run "env-C240-12x20"      INPUT-R12 12 20 4 NONE 2

echo "=== 20 cores, one socket against four ==="
ice_run "sock-C20-1x20-packed" INPUT    1 20  1 "$(packed 20)" 2
ice_run "sock-C20-1x20-spread" INPUT    1 20  1 "$(spread 20)" 2

echo "=== many blocks per rank: the R24 cut on one rank ==="
ice_run "blk-C80-1x80-R24"    INPUT-R24 1 80  1 "$(spread 80)" 2

echo "ICE-REDUCED DONE"
