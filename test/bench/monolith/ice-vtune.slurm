#!/usr/bin/env bash
#SBATCH --job-name=ice-vtune
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The profiler leg: the two single-rank configurations the 2026-09-30 wave left
# open -- one rank of 80 threads on the 192^3 block (70.9 % core-cycle efficiency
# against 82.1 % for 4 x 20) and the 24-block cut on one rank (per-cell phases
# 2-8x slower than with 4 blocks) -- plus one socket of 20 threads as the healthy
# reference. VTune wraps the binary through the harness's ICE_WRAP hook.
# Collections: hotspots (software sampling: the compute nodes' perf_event_paranoid
# may refuse hardware sampling), threading (spin and overhead per region), and
# memory-access, which needs hardware counters and is allowed to fail. Each
# result is reduced to text reports in $W/vt/<tag>/ (summary, hotspots by
# function, hotspots by source line of Mod_Fluxes).
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/vtune
source $ROOT/scripts/_ice_common.sh
ice_setup
source /opt/share/libs/intel/vtune/latest/env/vars.sh
vtune --version | head -1
mkdir -p $W/vt

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

# profile <tag> <set> <threads> <collect> [knobs...]
profile () {
  local tag=$1 set=$2 T=$3 col=$4; shift 4
  local r=$W/vt/$tag-$col
  rm -rf $r
  ICE_WRAP="vtune -collect $col $* -r $r -quiet --" \
  ice_run "vt-$tag-$col" $set 1 $T 1 "$(spread $T)" 1 || echo "VTUNE $tag $col: run failed (rc=$?)"
  if [ -d $r ]; then
    vtune -report summary  -r $r -format text -report-output $r.summary.txt  2>&1 | tail -2
    vtune -report hotspots -r $r -format text -group-by function -column "CPU Time" -report-output $r.functions.txt 2>&1 | tail -2
    vtune -report hotspots -r $r -format text -group-by source-line -source-object function=compute_flux -report-output $r.flux-lines.txt 2>&1 | tail -2
    [ $col = threading ] && vtune -report hotspots -r $r -format text -group-by function -column "Spin Time,Overhead Time,Effective Time" -report-output $r.spin.txt 2>&1 | tail -2
    echo "VTUNE $tag $col: reports in $r.*.txt"
  fi
}

echo "=== plain timing of the three legs, no profiler ==="
ice_run "plain-1x80"      INPUT     1 80 1 "$(spread 80)" 1
ice_run "plain-1x20"      INPUT     1 20 1 "$(spread 20)" 1
ice_run "plain-1x80-R24"  INPUT-R24 1 80 1 "$(spread 80)" 1

for col in hotspots threading; do
  echo "=== $col ==="
  profile 1x80      INPUT     80 $col -knob sampling-mode=sw
  profile 1x20      INPUT     20 $col -knob sampling-mode=sw
  profile 1x80-R24  INPUT-R24 80 $col -knob sampling-mode=sw
done

echo "=== memory-access (hardware counters: may be refused on this node) ==="
profile 1x80      INPUT     80 memory-access -knob analyze-mem-objects=false
profile 1x80-R24  INPUT-R24 80 memory-access -knob analyze-mem-objects=false

echo "ICE-VTUNE DONE"
