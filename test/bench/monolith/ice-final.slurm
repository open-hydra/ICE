#!/bin/bash
#SBATCH --job-name=ice-final
#SBATCH --partition=free
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Strong-scaling envelope for ICE, 1 -> 240 cores, plus the mode and layout
# sweeps. Marco Grossi's matrix of 2026-09-25 (jobs 302763..302995), unchanged;
# only the paths and the partition differ (README.md).
#
# The whole curve is taken in ONE job on ONE node set, because monolith's
# all-core clock varies ~4% between nodes and every efficiency is referred to
# this job's own 1-core anchor. Three nodes are held throughout even while the
# single-node points run: that wastes node-hours, but it is the only way the
# 240-core point shares an anchor with the 4-core one.
#
# No ballast. The MOSE campaign loaded the idle cores to hold the clock down;
# hydra-MF found that produced parallel efficiencies above 100% because the
# spinner is more power-hungry than the solver and pushed the anchor BELOW the
# full-node clock. ICE reports core-cycles per iteration instead, which is
# invariant under frequency scaling, so the clock divides out of C(1)/C(N)
# without perturbing the machine at all. Seconds are reported too, as what a
# user actually gets on an idle node.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/final
source $ROOT/scripts/_ice_common.sh
ice_setup

# spread(N): N cores taken evenly from all four sockets, matching how Intel MPI
# places ranks. Socket coverage is the first-order layout variable on this
# machine -- see the MOSE launch guide.
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

# packed(N): N cores taken in order, so a small N sits inside one socket. Used
# only for the socket-coverage comparison.
packed () { python3 -c "print(','.join(str(i) for i in range($1)))"; }

echo "=== 1. envelope: hybrid, 4 ranks per node, 1 -> 240 cores ==="
#           cores ranks thr ppn  set        reps
for spec in "1     1     1   1   INPUT      2" \
            "4     4     1   4   INPUT-R4   2" \
            "8     4     2   4   INPUT-R4   2" \
            "16    4     4   4   INPUT-R4   2" \
            "20    4     5   4   INPUT-R4   2" \
            "40    4    10   4   INPUT-R4   2" \
            "80    4    20   4   INPUT-R4   3" ; do
  set -- $spec; C=$1; R=$2; T=$3; P=$4; S=$5; REP=$6
  ice_run "env-C${C}-${R}x${T}" "$S" "$R" "$T" "$P" "$(spread $C)" "$REP"
done

echo "=== 2. envelope: multi-node, same per-node layout (4 ranks x 20 threads) ==="
ice_run "env-C160-8x20"  INPUT-R8   8 20 4 NONE 3
ice_run "env-C240-12x20" INPUT-R12 12 20 4 NONE 3

echo "=== 3. the three modes at matched core counts ==="
# Pure MPI stops where the node's memory does: every rank holds the whole
# domain, so 40 ranks would need ~220 GB of a 190 GB node.
# 20 ranks is absent on purpose: 192 = 2^6 x 3 has no factor of 5, so 192^3
# admits no balanced 20-way box decomposition -- MDB's best is 80% of ideal.
# 16 and 24 bracket it.
for spec in "2   INPUT-R2   2" \
            "4   INPUT-R4   2" \
            "8   INPUT-R8   2" \
            "16  INPUT-R16  2" \
            "24  INPUT-R24  2" ; do
  set -- $spec; C=$1; S=$2; REP=$3
  ice_run "mpi-C${C}-${C}x1" "$S" "$C" 1 "$C" "$(spread $C)" "$REP"
done
for C in 2 4 8 16 20 40 80 ; do
  ice_run "omp-C${C}-1x${C}" INPUT 1 "$C" 1 "$(spread $C)" 2
done

echo "=== 4. layout sweep at 80 cores, one node, full ==="
ice_run "lay-C80-4x20"  INPUT-R4   4 20  4 NONE 2
ice_run "lay-C80-8x10"  INPUT-R8   8 10  8 NONE 2
ice_run "lay-C80-16x5"  INPUT-R16 16  5 16 NONE 2
# 20x4 and 40x2 would complete the sweep, but 192^3 has no balanced 20- or
# 40-way box decomposition (no factor of 5), and 80 ranks would need ~500 GB.

echo "=== 5. socket coverage at 20 cores: one socket against all four ==="
ice_run "sock-C20-1x20-packed" INPUT     1 20 1 "$(packed 20)" 2
ice_run "sock-C20-1x20-spread" INPUT     1 20 1 "$(spread 20)" 2
ice_run "sock-C20-4x5-spread"  INPUT-R4  4  5 4 "$(spread 20)" 2

echo "=== 6. many blocks per rank: the R24 cut on one rank and on four ==="
# The campaign never ran more than one block per rank. ATLAS meshes and hydra
# cases do, and every per-block barrier of the solver (about 12 per stage) is
# paid once per block. The difference between these points and env-C80-4x20 /
# omp-C80-1x80 at equal core count is that cost plus the cut's ghost cells.
ice_run "blk-C80-1x80-R24" INPUT-R24  1 80  1 "$(spread 80)" 2
ice_run "blk-C80-4x20-R24" INPUT-R24  4 20  4 NONE 2

echo "=== 7. shared anchor, carried by every job in this campaign ==="
ice_run "anchor-4x20" INPUT-R4 4 20 4 NONE 2

echo "ICE-FINAL DONE"
