#!/usr/bin/env bash
#SBATCH --job-name=ice-transient2
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Second transient leg. ice-transient (303578) showed the start-up crawl on
# every binary that runs ANY loop through the OpenMP runtime's dispatcher
# (SCHEDULE(RUNTIME) or DYNAMIC, on the cell loops or on the tile loops alike)
# and on none without one; in the slow windows the whole process crawls, the
# memory-streaming loops 30-40x, flux 1.5-9x, while /usr/bin/time counts
# nothing (minor faults 47-96 k, involuntary switches 4-11 k per run). The
# kernel's NUMA balancer is on and its first scan pass over the 7 GB footprint
# lasts about as long as the crawl (~45-60 s of wall time, whatever the
# per-iteration cost), yet the old binary at the 24-block cut takes 140 k
# hinting faults and stays flat: balancer activity alone is not the cost.
#
# This job: V3 = C5d's shallow tiling with no dispatcher anywhere (the design
# control: is a dispatcher-free kernel flat, and what does it cost at steady
# state against the dynamic tiles?); the instrument binary with its pages
# interleaved over the four nodes from the start (numactl --interleave=all:
# if the crawl needs the serial first touch's node-0 concentration, this is
# flat), with KMP_BLOCKTIME=infinite (no thread ever sleeps: rules out the
# futex path), and C5d interleaved (does the gain survive). Every leg samples
# /proc/<pid>/task/<tid>/sched of four threads and numastat -p every 5 s: a
# non-zero numa_group_id says the kernel grouped the threads on a shared page
# and herds their memory toward one node.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/transient2
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

echo "=== kernel knobs on $(hostname) ==="
for f in /proc/sys/kernel/numa_balancing /proc/sys/kernel/numa_balancing_scan_delay_ms \
         /proc/sys/kernel/numa_balancing_scan_period_min_ms /proc/sys/kernel/numa_balancing_scan_period_max_ms \
         /proc/sys/kernel/numa_balancing_scan_size_mb /sys/kernel/mm/transparent_hugepage/enabled \
         /sys/kernel/mm/transparent_hugepage/defrag; do
  [ -r $f ] && echo "$f: $(cat $f)"
done
uname -r; which numactl || echo "no numactl"

VM='^(numa_hint_faults|numa_hint_faults_local|numa_pages_migrated|numa_pte_updates|numa_huge_pte_updates|pgmigrate_success|pgmigrate_fail|thp_fault_alloc|thp_fault_fallback|thp_split_page|thp_deferred_split_page|pgfault) '
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

# numa_sampler <tag>: every 5 s, for the rank process (basename ICE-mpi), the
# per-node memory (numastat -p) and the NUMA-balancing state of four threads.
numa_sampler () {
  local out=$W/logs/$1-numa.log; : > $out
  while :; do
    local pid=$(pgrep -x ICE-mpi | head -1)
    if [ -n "$pid" ] && [ -d /proc/$pid ]; then
      echo "T=$SECONDS pid=$pid" >> $out
      numastat -p $pid 2>/dev/null | grep -E '^(Total|Private|Heap|Stack)' >> $out
      for t in $(ls /proc/$pid/task 2>/dev/null | sort -n | awk 'NR==1||NR==21||NR==41||NR==61'); do
        echo "task $t: $(grep -E 'numa_pages_migrated|numa_preferred_nid|total_numa_faults|numa_group_id|^numa_faults' /proc/$pid/task/$t/sched 2>/dev/null | tr -s ' ' | tr '\n' ';')" >> $out
      done
    fi
    sleep 5
  done
}

# leg <tag> <binary> <set> [wrap] [env...]: one run of 60 iterations at 1 x 80.
leg () {
  local tag=$1 bin=$2 set=$3 wrap=${4:-}; shift 3; [ $# -gt 0 ] && shift
  [ -x "$bin" ] || { echo "LEG $tag: binary $bin missing"; return; }
  export ICE_PERF_BIN=$bin
  source $ROOT/scripts/_ice_common.sh
  ice_setup
  export ICE_WRAP="$wrap /usr/bin/time -v"
  local e; for e in "$@"; do export "$e"; done
  echo "LEG tag=$tag bin=$bin set=$set wrap='$wrap' env='$*'"
  numa_sampler "$tag" & local spid=$!
  local before=$(vm_snap)
  ice_run "$tag" $set 1 80 1 "$(spread 80)" 1
  local after=$(vm_snap)
  kill $spid 2>/dev/null; wait $spid 2>/dev/null
  for e in "$@"; do unset "${e%%=*}"; done
  vm_delta "$before" "$after"
  local log=$W/logs/${tag}-r1.log
  echo "WINDOWS $tag: $(grep -o 'wall/iter *[0-9.E+-]*' $log | awk '{printf " %.3f", $2}')  avg $(grep -o 'Solver per iteration *[0-9.E+-]*' $log | awk '{print $4}')"
  grep -E 'Minor \(reclaiming|Voluntary context|Involuntary context' $log | sed -E 's/^\s*/  TIME /'
  echo "NUMAGROUP $tag: $(grep -o 'numa_group_id=[0-9]*' $W/logs/$tag-numa.log | sort | uniq -c | tr '\n' ' ')"
  echo "PREFNID $tag: $(grep -o 'numa_preferred_nid : -\?[0-9]*' $W/logs/$tag-numa.log | sort | uniq -c | tr -s ' ' | tr '\n' ' ')"
}

OLD=$ROOT/final/bin.303485/ICE-mpi
INST=$ROOT/sched/bin.303564/ICE-mpi
C5D=$ROOT/bin/ICE-mpi
V3=$ROOT/src-v3/bin-hyb/ICE
for b in $OLD $INST $C5D $V3; do echo "BIN $b $(md5sum $b 2>/dev/null | cut -c1-12)"; done

echo "=== 1 x 80, one 192^3 block, 60 iterations ==="
leg t2-old-1x80        $OLD  INPUT
leg t2-inst-1x80       $INST INPUT
leg t2-v3-1x80         $V3   INPUT
leg t2-inst-interleave $INST INPUT "numactl --interleave=all"
leg t2-inst-nosleep    $INST INPUT "" KMP_BLOCKTIME=infinite
leg t2-c5d-interleave  $C5D  INPUT "numactl --interleave=all"

echo "=== the 24-block cut ==="
leg t2-v3-R24          $V3   INPUT-R24
leg t2-c5d-interleave-R24 $C5D INPUT-R24 "numactl --interleave=all"

echo "ICE-TRANSIENT2 DONE"
