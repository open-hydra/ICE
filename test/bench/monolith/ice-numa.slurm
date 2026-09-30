#!/usr/bin/env bash
#SBATCH --job-name=ice-numa
#SBATCH --partition=free
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --output=/data_fast/gpassarani/ice-perf/logs/%x-%j.out
#
# The single-rank legs after the OpenMP wave: one rank of 40 or 80 threads over
# two or four sockets, with the pages where the serial first touch left them
# (plain) and interleaved over the sockets (inter). On the baseline binary
# (303351) interleaving lost 7-9 %; the gather kernel is faster per face and
# therefore more exposed to where the pages are (303387/303402: 1x80 slower than
# the baseline while 4x20 is faster). The answer decides the first-touch phase
# (C3). No anchor: the question is the ratio inside this job.
set -u
ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
W=$ROOT/numa
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

for T in 80 40; do
  echo "=== 1 x $T: pages where the serial first touch left them, then interleaved ==="
  ice_run "plain-1x${T}"  INPUT 1 $T 1 "$(spread $T)" 2
  ICE_WRAP="numactl --interleave=all" \
  ice_run "inter-1x${T}"  INPUT 1 $T 1 "$(spread $T)" 2
done

echo "=== the four-rank layout on the same node, for the ratio ==="
ice_run "env-C80-4x20"    INPUT-R4 4 20 4 NONE 2

echo "ICE-NUMA DONE"
