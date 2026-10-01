#!/bin/bash
#SBATCH --job-name=ice-smoke
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# Pre-flight for the ICE campaign: does the binary run here, does the pinning
# take, how long is set-up, and how much memory does a rank really hold?
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/smoke
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

lscpu | grep -E "^Socket|^Core|^CPU\(s\)|Model name|Thread" || true
free -g | head -2
cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null | sed 's/^/perf_event_paranoid /'

ice_run "smoke-1x4"   INPUT     1  4 1 "$(spread 4)"  1
ice_run "smoke-4x1"   INPUT-R4  4  1 4 "$(spread 4)"  1
ice_run "smoke-4x20"  INPUT-R4  4 20 4 NONE           1
ice_run "smoke-16x1"  INPUT-R16 16 1 16 "$(spread 16)" 1

echo "ICE-SMOKE DONE"
