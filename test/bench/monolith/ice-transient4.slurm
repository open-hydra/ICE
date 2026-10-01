#!/usr/bin/env bash
#SBATCH --job-name=ice-transient4
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:10:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The magnitude question. ice-transient2 (303579) showed the crawl happening
# BEFORE the balancer moves anything (node 0 holds 4.4 GB unchanged through
# the first two slow windows; the migration of 1.7 GB comes in the third) and
# vanishing when the balancer is off (numactl --interleave=all). Sixty remote
# threads streaming one node's memory should cost 2-3x, not the 6-10x
# measured, so two readings of the same data compete: (a) the kernel's
# hinting/scan pass on an 80-thread THP footprint is what crawls, the
# placement only its occasion; (b) node 0's memory system saturates and the
# crawl IS the placement. numactl --membind=0 decides: the whole footprint on
# node 0 for the entire run, the balancer off (explicit policy). Flat at
# ~0.4 s -> (a), remote access is a 10-15 % tax and removing the trigger is
# enough; 2 s per iteration throughout -> (b), first-touch placement is
# mandatory for any single-rank number to mean anything. The old binary runs
# the same leg for the tax it pays.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/transient4
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
OLDR=$ROOT/src-old/bin-hyb/ICE        # b151bea's source rebuilt 2026-10-01 with today's build-cfg.sh
INST=$ROOT/sched/bin.303564/ICE-mpi
for b in $OLD $OLDR $INST; do echo "BIN $b $(md5sum $b 2>/dev/null | cut -c1-12)"; done

# ice-transient3 (303580) found the master pinned to CPU 0 through set-up in EVERY binary, and
# both single-factor variants crawling (V5: old kernel + the Wrap_Setup block; V6: shallow
# static kernel + old Wrap_Setup): two unrelated edits cannot both be "the" trigger. If the
# same source as the flat snapshot, rebuilt today, crawls, the trigger is the build, not
# the code.
echo "=== 1 x 80: the flat snapshot against its own source rebuilt today ==="
leg t4-old-1x80    $OLD  INPUT
leg t4-oldr-1x80   $OLDR INPUT
leg t4-oldr2-1x80  $OLDR INPUT

echo "=== 1 x 80, one 192^3 block, 60 iterations, the whole footprint bound to node 0 ==="
leg t4-inst-membind0 $INST INPUT "numactl --membind=0"
leg t4-old-membind0  $OLD  INPUT "numactl --membind=0"

echo "ICE-TRANSIENT4 DONE"
