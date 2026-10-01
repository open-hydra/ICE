#!/usr/bin/env bash
#SBATCH --job-name=ice-cuts
#SBATCH --partition=free
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The decomposition leg: the greedy cuts of the campaign (INPUT-R<n>, ATLAS MDB
# with a hand-swept min-cells) against the halo-objective cuts of MDB's search
# (INPUT-R<n>-halo, ATLAS branch feat/mdb-optimal-decomposition), at equal rank
# and thread counts, so that only the shape of the blocks differs. Same protocol
# as the other jobs: one job, one node set, the job's own 1-core anchor, minimum
# over repetitions, core-cycles beside seconds. The ghost share of each set comes
# from decompositions.tsv (analyze.py's third argument).
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/cuts
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

echo "=== anchor ==="
ice_run "env-C1-1x1"            INPUT          1  1  1 "$(spread 1)"  2

echo "=== R12 on three nodes: the halo at 240 cores ==="
ice_run "cut-C240-12x20-R12"    INPUT-R12     12 20  4 NONE 2
ice_run "cut-C240-12x20-R12h"   INPUT-R12-halo 12 20 4 NONE 2

echo "=== R16: pure MPI on one node, and 16 x 5 ==="
ice_run "cut-C16-16x1-R16"      INPUT-R16     16  1 16 "$(spread 16)" 2
ice_run "cut-C16-16x1-R16h"     INPUT-R16-halo 16 1 16 "$(spread 16)" 2
ice_run "cut-C80-16x5-R16"      INPUT-R16     16  5 16 NONE 2
ice_run "cut-C80-16x5-R16h"     INPUT-R16-halo 16 5 16 NONE 2

echo "=== R20: the count the greedy cut could not balance ==="
ice_run "cut-C20-20x1-R20h"     INPUT-R20-halo 20 1 20 "$(spread 20)" 2

echo "=== R24: pure MPI, 4 x 20 with six blocks per rank, and all 24 on one rank ==="
ice_run "cut-C24-24x1-R24"      INPUT-R24     24  1 24 "$(spread 24)" 2
ice_run "cut-C24-24x1-R24h"     INPUT-R24-halo 24 1 24 "$(spread 24)" 2
ice_run "cut-C80-4x20-R24"      INPUT-R24      4 20  4 NONE 2
ice_run "cut-C80-4x20-R24h"     INPUT-R24-halo 4 20  4 NONE 2
ice_run "cut-C80-1x80-R24"      INPUT-R24      1 80  1 "$(spread 80)" 2
ice_run "cut-C80-1x80-R24h"     INPUT-R24-halo 1 80  1 "$(spread 80)" 2

echo "ICE-CUTS DONE"
