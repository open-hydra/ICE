#!/bin/bash
#SBATCH --job-name=ice-attr
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=05:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Attribution of the one-group speed-up (plan openmp-loose-ends, L1). Job 303832
# read, with ONE thread group, 0.501 s on the 24-block cut at 1 x 80 against 0.550
# for the phase-2 binary (935af6b), and 0.295 against 0.305 on the single block --
# best of two, while two runs of one layout differed by 9 % in the same job.
# Candidates: the merge of main (037e94f, set-up only), the per-cell phases by
# group (e69ebcc), ghost/bound bodies in procedures (fddd1cb), the flux phases in
# tile procedures (055f6d4), the tip (f50efc2) -- or noise.
#
# Every binary built by build-attr-monolith.sh from one clone; ICE_THREAD_GROUPS=1
# everywhere (the older binaries ignore it); interleaved order, so drift of the
# node hits every binary alike: FIVE passes of (each binary: R24 1x80, single
# block 1x80). attr-summary.py reads the IOK lines (steady_iter, ICE Detail).
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/attr/run
B=$ROOT/attr/bin
COMMITS=${ICE_ATTR_COMMITS:-935af6b 037e94f e69ebcc fddd1cb 055f6d4 f50efc2}
PASSES=${ICE_ATTR_PASSES:-5}
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
for c in $COMMITS; do
  [ -x $B/ICE-mpi.$c ] || { echo "binary $B/ICE-mpi.$c missing"; exit 2; }
  echo "BIN $c $(md5sum $B/ICE-mpi.$c | cut -c1-12) $(head -1 $B/ICE-mpi.$c.commit)"
done
export ICE_THREAD_GROUPS=1
for k in $(seq 1 $PASSES); do
  echo "=== pass $k ==="
  for c in $COMMITS; do
    export ICE_PERF_BIN=$B/ICE-mpi.$c
    source $ROOT/scripts/_ice_common.sh
    ice_setup
    ice_run "a-$c-R24-p$k" INPUT-R24 1 80 1 "$(spread 80)" 1
    ice_run "a-$c-R1-p$k"  INPUT     1 80 1 "$(spread 80)" 1
  done
done
echo "ICE-ATTR DONE"
