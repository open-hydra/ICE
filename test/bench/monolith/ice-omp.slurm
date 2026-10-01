#!/bin/bash
#SBATCH --job-name=ice-omp
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Pure-OpenMP placement sweep: one socket against four, and page interleaving.
# Marco Grossi's job 302996 (2026-09-28), unchanged but for the paths.
#
# Its result on the campaign binary (AUDIT.md, section 3f): spreading twenty
# threads over four sockets was FASTER than packing them on the socket that
# holds their pages (95.7 % against 90.6 % cycle efficiency), and --interleave
# made 1x40 and 1x80 slower. The kernel is dominated by indirect calls, not by
# memory traffic, so page placement does not bind. This job is re-run after the
# flux kernel and the dispatch are fixed: once the kernel is memory-bound the
# same three series decide whether a parallel first touch is worth doing.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/omp
source $ROOT/scripts/_ice_common.sh
ice_setup

# socket0(N): the first N cores of socket 0. Confines everything -- threads,
# barriers and memory -- to one NUMA domain.
sock0 () { python3 - "$1" <<'PY'
import subprocess, collections, sys
N=int(sys.argv[1])
out=subprocess.run(['lscpu','-p=CPU,SOCKET'],stdout=subprocess.PIPE).stdout.decode()
sock=collections.OrderedDict()
for line in out.splitlines():
    if line.startswith('#') or not line.strip(): continue
    c,s=line.split(',')[:2]; sock.setdefault(int(s),[]).append(int(c))
for s in sock: sock[s].sort()
print(','.join(map(str,sock[list(sock)[0]][:N])))
PY
}

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

echo "=== A. pure OpenMP INSIDE one socket: NUMA-free, barrier-local ==="
# The anchor is core 0 of socket 0, so every point in this series shares both
# the socket and the reference. Scaling here is ICE's OpenMP quality with page
# placement held innocent by construction.
for T in 1 2 4 5 10 20 ; do
  ice_run "s0-1x${T}" INPUT 1 $T 1 "$(sock0 $T)" 2
done

echo "=== B. the same thread counts SPREAD over four sockets ==="
# Identical work, identical thread count, identical binary; only the pages are
# now remote for three quarters of the threads. The difference between this
# series and A at the same T is the cost of page placement alone.
for T in 4 8 20 40 ; do
  ice_run "sp-1x${T}" INPUT 1 $T 1 "$(spread $T)" 2
done

echo "=== C. does interleaving the pages recover it? ==="
# --interleave is not the fix -- it makes every page equally remote instead of
# local -- but it removes the single-controller bottleneck, so a large gain
# here is positive evidence that first-touch placement is what costs the 1x80
# point its efficiency.
for T in 40 80 ; do
  ice_run "plain-1x${T}"  INPUT 1 $T 1 "$(spread $T)" 2
  ICE_WRAP="numactl --interleave=all" \
  ice_run "inter-1x${T}"  INPUT 1 $T 1 "$(spread $T)" 2
done

echo "=== D. hybrid controls at 80 cores, for the same reference ==="
# 4x20 puts one rank on each socket and 16x5 puts four; both are NUMA-clean by
# construction because a rank first-touches its own pages on its own socket.
# If they still differ, page placement is not the whole story and the residue
# is threads-per-rank.
ice_run "hyb-4x20"  INPUT-R4  4 20 4 NONE 2
ice_run "hyb-16x5"  INPUT-R16 16 5 16 NONE 2

echo "ICE-OMP DONE"
