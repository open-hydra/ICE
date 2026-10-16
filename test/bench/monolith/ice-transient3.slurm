#!/usr/bin/env bash
#SBATCH --job-name=ice-transient3
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Third transient leg. ice-transient2 (303579) placed the crawl: the binaries
# that crawl start the loop with their serially first-touched memory ALL on
# node 0 (numastat during set-up: 2952 MB on node 0, 7/2/19 MB elsewhere),
# the old binary with it already spread (2234 / 475 / 348 / 376 MB at the same
# point of set-up); the crawl lasts until the kernel's balancer has moved
# ~1.7 GB off node 0, after which both reach the same distribution and speed.
# Under numactl --interleave=all (balancer off, pages spread) every binary is
# flat. The set-up code is identical, so what differs is WHERE THE MASTER
# THREAD RUNS during read_ic/allocate_data: stationary on node 0 or wandering.
# This job traces it: the master's Cpus_allowed_list and current CPU every
# second through set-up, KMP_AFFINITY=verbose on the first two legs (the
# runtime prints each binding as it happens, in order with ICE's own output),
# and the two single-factor variants V5 (old flux kernel + the Wrap_Setup
# block with omp_set_schedule) and V6 (shallow static kernel, old Wrap_Setup),
# plus V4 (shallow kernel, SCHEDULE(STATIC,1) round-robin tiles, no dispatcher).
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/transient3
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
master_trace () {   # every second: the master thread's allowed CPUs, current CPU and node, and the thread count
  local out=$W/logs/$1-master.log; : > $out
  while :; do
    local pid=$(pgrep -x ICE-mpi | head -1)
    if [ -n "$pid" ] && [ -d /proc/$pid ]; then
      local cpu=$(awk '{print $39}' /proc/$pid/task/$pid/stat 2>/dev/null)
      local node=$(cat /sys/devices/system/cpu/cpu${cpu:-0}/topology/physical_package_id 2>/dev/null)
      echo "T=$SECONDS cpu=$cpu pkg=$node allowed=$(grep Cpus_allowed_list /proc/$pid/task/$pid/status | awk '{print $2}') threads=$(ls /proc/$pid/task | wc -l) rss_mb=$(awk '/VmRSS/{printf "%d", $2/1024}' /proc/$pid/status)" >> $out
    fi
    sleep 1
  done
}
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
  master_trace "$tag" & local mpid2=$!
  local before=$(vm_snap)
  ice_run "$tag" $set 1 80 1 "$(spread 80)" 1
  local after=$(vm_snap)
  kill $spid $mpid2 2>/dev/null; wait $spid $mpid2 2>/dev/null
  for e in "$@"; do unset "${e%%=*}"; done
  vm_delta "$before" "$after"
  local log=$W/logs/${tag}-r1.log
  echo "WINDOWS $tag: $(grep -o 'wall/iter *[0-9.E+-]*' $log | awk '{printf " %.3f", $2}')  avg $(grep -o 'Solver per iteration *[0-9.E+-]*' $log | awk '{print $4}')"
  grep -E 'Minor \(reclaiming|Voluntary context|Involuntary context' $log | sed -E 's/^\s*/  TIME /'
  echo "MASTER $tag: $(awk '{printf "%s/%s ", $3, $2}' $W/logs/$tag-master.log | cut -c1-400)"
  echo "ALLOWED $tag: $(awk '{print $4}' $W/logs/$tag-master.log | uniq -c | tr -s ' ' | tr '\n' ' ')"
  grep -E 'OMP: Info #(171|247|242)|bound to OS proc set' $log | head -3 | sed 's/^/  KMP /'
  echo "NUMAGROUP $tag: $(grep -o 'numa_group_id=[0-9]*' $W/logs/$tag-numa.log | sort | uniq -c | tr '\n' ' ')"
  echo "PREFNID $tag: $(grep -o 'numa_preferred_nid : -\?[0-9]*' $W/logs/$tag-numa.log | sort | uniq -c | tr -s ' ' | tr '\n' ' ')"
}

OLD=$ROOT/final/bin.303485/ICE-mpi
INST=$ROOT/sched/bin.303564/ICE-mpi
V4=$ROOT/src-v4/bin-hyb/ICE
V5=$ROOT/src-v5/bin-hyb/ICE
V6=$ROOT/src-v6/bin-hyb/ICE
for b in $OLD $INST $V4 $V5 $V6; do echo "BIN $b $(md5sum $b 2>/dev/null | cut -c1-12)"; done

echo "=== 1 x 80, one 192^3 block, 60 iterations ==="
leg t3-old-1x80  $OLD  INPUT "" KMP_AFFINITY=verbose
leg t3-inst-1x80 $INST INPUT "" KMP_AFFINITY=verbose
leg t3-v5-1x80   $V5   INPUT
leg t3-v6-1x80   $V6   INPUT
leg t3-v4-1x80   $V4   INPUT

echo "ICE-TRANSIENT3 DONE"
