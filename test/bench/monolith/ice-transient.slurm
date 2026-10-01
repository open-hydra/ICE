#!/usr/bin/env bash
#SBATCH --job-name=ice-transient
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The start-up transient of the SCHEDULE(RUNTIME) binaries (found 2026-10-01).
# Every run of the instrument binary 3ad9b2b3 at 1 x 80 (jobs 303564, 303570)
# starts at 2.1-2.5 s per iteration and reaches its steady 0.35-0.45 s only
# after 20-30 iterations; the binaries before it run flat from the first window
# (0.42 -> 0.36). The harness's wall_iter is the last window, so the transient
# went unseen. Its slow windows are the per-cell collapse(3) loops (zero 333 ms
# against 8, copy 272 against 8, source 634 against 41) while the tile loops,
# through the same dispatcher, cost 1.35x -- a page-migration time scale.
#
# This job bisects the trigger on one node with 60 iterations per run (six
# timer windows): the pre-instrument binary (snapshot of job 303485), the C5d
# binary (cell loops schedule(runtime) set static by ICE, tile loops
# DYNAMIC,1, shallow tiling), V1 = C5d with the cell loops back to compile-time
# static, V2 = C5d with the tile loops back to STATIC and the cell loops still
# runtime, and the instrument binary itself (snapshot of job 303564), the one
# known to show the transient: without it a flat node is indistinguishable from
# a cured trigger. The kernel's NUMA-balancing and THP counters are read around each
# run (/proc/vmstat deltas) and /usr/bin/time -v wraps the binary for the
# page-fault and context-switch counts. The same four at the 24-block cut.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/transient
mkdir -p $W
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
         /sys/kernel/mm/transparent_hugepage/defrag /sys/kernel/mm/transparent_hugepage/khugepaged/pages_to_scan \
         /sys/kernel/mm/transparent_hugepage/khugepaged/scan_sleep_millisecs; do
  [ -r $f ] && echo "$f: $(cat $f)"
done
uname -r; numactl -H | grep -E 'available|node [0-9] cpus' | head -5

VM='^(numa_hint_faults|numa_hint_faults_local|numa_pages_migrated|numa_pte_updates|numa_huge_pte_updates|pgmigrate_success|pgmigrate_fail|thp_fault_alloc|thp_fault_fallback|thp_collapse_alloc|thp_split_page|thp_deferred_split_page|pgfault|compact_stall) '
vm_snap () { grep -E "$VM" /proc/vmstat | tr '\n' ' '; echo; }
vm_delta () {   # vm_delta <before-line> <after-line>
  python3 - "$1" "$2" <<'PY'
import sys
def parse(s):
    t=s.split(); return {t[i]:int(t[i+1]) for i in range(0,len(t)-1,2)}
a=parse(sys.argv[1]); c=parse(sys.argv[2])
print('VMDELTA ' + ' '.join('%s=%d' % (k, c[k]-a[k]) for k in a if k in c))
PY
}

# leg <tag> <binary> <set>: one run of 60 iterations at 1 x 80, counters around it.
leg () {
  local tag=$1 bin=$2 set=$3
  [ -x "$bin" ] || { echo "LEG $tag: binary $bin missing"; return; }
  export ICE_PERF_BIN=$bin
  source $ROOT/scripts/_ice_common.sh
  ice_setup
  if [ -x /usr/bin/time ]; then export ICE_WRAP="/usr/bin/time -v"; else unset ICE_WRAP; echo "  (no /usr/bin/time on this node: no fault counts)"; fi
  echo "LEG tag=$tag bin=$bin set=$set"
  local before=$(vm_snap)
  ice_run "$tag" $set 1 80 1 "$(spread 80)" 1
  local after=$(vm_snap)
  vm_delta "$before" "$after"
  local log=$W/logs/${tag}-r1.log
  echo "WINDOWS $tag: $(grep -o 'wall/iter *[0-9.E+-]*' $log | awk '{printf " %.3f", $2}')  avg $(grep -o 'Solver per iteration *[0-9.E+-]*' $log | awk '{print $4}')"
  grep -E 'Minor \(reclaiming|Major \(requiring|Voluntary context|Involuntary context|Maximum resident|Elapsed' $log | sed 's/^/  TIME /'
}

OLD=$ROOT/final/bin.303485/ICE-mpi
INST=$ROOT/sched/bin.303564/ICE-mpi   # the instrument binary 3ad9b2b3, the one known to show the transient
C5D=$ROOT/bin/ICE-mpi
V1=$ROOT/src-v1/bin-hyb/ICE
V2=$ROOT/src-v2/bin-hyb/ICE
for b in $OLD $INST $C5D $V1 $V2; do echo "BIN $b $(md5sum $b 2>/dev/null | cut -c1-12)"; done

echo "=== 1 x 80, one 192^3 block, 60 iterations ==="
leg tr-old-1x80 $OLD INPUT
leg tr-inst-1x80 $INST INPUT
leg tr-c5d-1x80 $C5D INPUT
leg tr-v1-1x80  $V1  INPUT
leg tr-v2-1x80  $V2  INPUT

echo "=== 1 x 80, the 24-block cut ==="
leg tr-old-R24 $OLD INPUT-R24
leg tr-inst-R24 $INST INPUT-R24
leg tr-c5d-R24 $C5D INPUT-R24
leg tr-v1-R24  $V1  INPUT-R24
leg tr-v2-R24  $V2  INPUT-R24

echo "ICE-TRANSIENT DONE"
