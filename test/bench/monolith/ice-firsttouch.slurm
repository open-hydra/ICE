#!/usr/bin/env bash
#SBATCH --job-name=ice-firsttouch
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The remedy leg. Jobs 303578-303745 traced the 1 x 80 start-up crawl to the
# block fields' serial first touch (Allocate_Block writes prim, prim_old, tau,
# beta, source, residual from the master, pinned on CPU 0: the whole block
# starts on node 0) met by the kernel's NUMA balancer, which needs ~45-60 s to
# move ~1.7 GB. 038824d (reverted 09-30 on a last-window metric, on a binary
# that did not crawl) writes the same values from a static parallel loop over
# the k-planes, so each plane starts on the socket of the thread that sweeps
# it. This job runs it on C5d (C5DFT) and on the old code (OLDFT) beside both
# without it, two repetitions each, single block and the 24-block cut.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/firsttouch
mkdir -p $W $W/logs
INI60=$W/bench60.ini
sed 's/^iter-threshold *= *40/iter-threshold = 60/' $ROOT/scripts/bench.ini > $INI60
grep -q 'iter-threshold = 60' $INI60 || { echo "bench60.ini: iter-threshold not rewritten"; exit 2; }
export ICE_PERF_INI=$INI60

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

echo "=== $(hostname): numa_balancing=$(cat /proc/sys/kernel/numa_balancing) THP=$(cat /sys/kernel/mm/transparent_hugepage/enabled) ==="
numactl -H | grep -E 'node [0-9] (size|free)'

VM='^(numa_hint_faults|numa_pages_migrated|numa_huge_pte_updates|thp_fault_alloc|pgfault) '
vm_snap () { grep -E "$VM" /proc/vmstat | tr '\n' ' '; echo; }
vm_delta () {
  python3 - "$1" "$2" <<'PY'
import sys
def parse(s):
    t=s.split(); return {t[i]:int(t[i+1]) for i in range(0,len(t)-1,2)}
a=parse(sys.argv[1]); c=parse(sys.argv[2])
print('VMDELTA ' + ' '.join('%s=%d' % (k, c[k]-a[k]) for k in a if k in c))
PY
}
mem_sampler () {   # numastat -p of the rank every 5 s
  local out=$W/logs/$1-numa.log; : > $out
  while :; do
    local pid=$(pgrep -x ICE-mpi | head -1)
    [ -n "$pid" ] && [ -d /proc/$pid ] && { echo "T=$SECONDS" >> $out; numastat -p $pid 2>/dev/null | grep -E '^Total' >> $out; }
    sleep 5
  done
}

leg () {   # leg <tag> <binary> <set> [wrap]
  local tag=$1 bin=$2 set=$3 wrap=${4:-}
  [ -x "$bin" ] || { echo "LEG $tag: binary $bin missing"; return; }
  export ICE_PERF_BIN=$bin
  source $ROOT/scripts/_ice_common.sh
  ice_setup
  export ICE_WRAP="$wrap /usr/bin/time -v"
  echo "LEG tag=$tag bin=$bin set=$set wrap='$wrap'"
  mem_sampler "$tag" & local spid=$!
  local before=$(vm_snap)
  ice_run "$tag" $set 1 80 1 "$(spread 80)" 1
  local after=$(vm_snap)
  kill $spid 2>/dev/null; wait $spid 2>/dev/null
  vm_delta "$before" "$after"
  local log=$W/logs/${tag}-r1.log
  echo "WINDOWS $tag: $(grep -o 'wall/iter *[0-9.E+-]*' $log | awk '{printf " %.3f", $2}')  avg $(grep -o 'Solver per iteration *[0-9.E+-]*' $log | awk '{print $4}')"
  echo "NODES $tag (MB on nodes 0-3, loop start / end): $(grep '^Total' $W/logs/$tag-numa.log | awk 'NR==8{printf "%s %s %s %s", $2,$3,$4,$5} END{printf " / %s %s %s %s", $2,$3,$4,$5}')"
}

OLD=$ROOT/final/bin.303485/ICE-mpi
C5D=$ROOT/bin/ICE-mpi
C5DFT=$ROOT/src-c5dft/bin-hyb/ICE
OLDFT=$ROOT/src-oldft/bin-hyb/ICE
for b in $OLD $C5D $C5DFT $OLDFT; do echo "BIN $b $(md5sum $b 2>/dev/null | cut -c1-12)"; done

echo "=== 1 x 80, one 192^3 block, 60 iterations ==="
for r in 1 2; do
  leg ft-c5dft-1x80-$r $C5DFT INPUT
  leg ft-c5d-1x80-$r   $C5D   INPUT
  leg ft-oldft-1x80-$r $OLDFT INPUT
  leg ft-old-1x80-$r   $OLD   INPUT
done
echo "=== 1 x 80, the 24-block cut ==="
leg ft-c5dft-R24 $C5DFT INPUT-R24
leg ft-c5d-R24   $C5D   INPUT-R24
leg ft-old-R24   $OLD   INPUT-R24

echo "ICE-FIRSTTOUCH DONE"
