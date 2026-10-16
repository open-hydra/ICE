# Shared harness for the ICE scaling jobs on monolith. Sourced, not executed.
#
# Origin: Marco Grossi's campaign of 2026-09-25..28 (jobs 302763..302995, the numbers
# of docs/user/parallel.md). Taken over on 2026-09-30: paths re-pointed to
# ICE_PERF_ROOT, every IOK field now comes from the BEST repetition, OMP_PROC_BIND
# is set, and the job log carries the ICE commit. See README.md beside this file.
#
# Conventions inherited from the MOSE and hydra-MF campaigns, each of which cost
# something to learn there:
#   * score the LAST timer window -- the first carries set-up transients;
#   * the decomposition belongs to the RANK count, never the core count;
#   * every point in a curve must be measured on the SAME node -- monolith's
#     all-core clock varies ~4% node to node;
#   * efficiency is referred to that node's own 1-core anchor, from the same job;
#   * repeat the decisive points and keep the minimum;
#   * report core-cycles as well as seconds. Seconds charge the solver for the
#     chip's frequency scaling (one active core turbos to ~3.8 GHz here while
#     eighty sit near ~2.6 GHz); total core-cycles per iteration are invariant
#     under perfect parallelism, so C(1)/C(N) is efficiency with the clock
#     already divided out. hydra-MF abandoned the MOSE ballast method for this
#     one after the ballast produced efficiencies above 100%.
#
# ICE-specific: every rank allocates the WHOLE domain (Wrap_Setup says so in as
# many words), so per-rank RSS does not fall as ranks are added and the node's
# memory, not its cores, is what caps pure MPI. RSS is therefore sampled for
# every point and published beside the timing.

ROOT=${ICE_PERF_ROOT:-/data_fast/gpassarani/ice-perf}
BIN=${ICE_PERF_BIN:-$ROOT/bin/ICE-mpi}
CASE=${ICE_PERF_CASE:-$ROOT/case}
INI=${ICE_PERF_INI:-$ROOT/scripts/bench.ini}

ice_setup () {
  ulimit -s unlimited
  export KMP_STACKSIZE=200M OMP_STACKSIZE=200M
  export I_MPI_PIN_RESPECT_CPUSET=1 I_MPI_PIN_DOMAIN=omp
  # OMP_PLACES=cores alone lets a thread migrate between the cores of its
  # rank's domain; close binding pins it. Never set in the 2026-09 campaign.
  export OMP_PROC_BIND=close
  mkdir -p $W/logs $W/work
  cd $W

  # SNAPSHOT the binary into the job's own directory and run that copy.
  #
  # A rebuild into $ROOT/bin while a job is running silently changes the
  # binary underneath it: mpirun execs the file afresh for every repetition, so
  # the points before the rebuild and the points after it come from different
  # code and the curve is a splice of two binaries that still looks like one.
  # That happened on 2026-09-28 and cost job 302988. A curve must be one
  # binary, and this is the only way to guarantee it from inside the job.
  # The copy keeps the BASENAME `ICE-mpi` and is isolated by DIRECTORY instead.
  # The RSS sampler below finds the ranks with `ps -C ICE-mpi`, which matches on
  # the command name, so a snapshot named `ICE-mpi.$SLURM_JOB_ID` matches
  # nothing and every rss_node_gb in the job comes back 0.00. That is what
  # happened to job 302995, whose memory column is empty for this reason.
  mkdir -p "$W/bin.$SLURM_JOB_ID"
  cp -f "$BIN" "$W/bin.$SLURM_JOB_ID/ICE-mpi" && BIN="$W/bin.$SLURM_JOB_ID/ICE-mpi"

  echo "JOB=$SLURM_JOB_NAME id=$SLURM_JOB_ID nodes=$SLURM_NODELIST partition=${SLURM_JOB_PARTITION:-?}"
  echo "BINARY md5=$(md5sum $BIN | cut -c1-12) $(ls -l $BIN | awk '{print $5}') bytes (snapshot: $BIN)"
  # The commit the binary was built from, written by sync-src-monolith.sh next
  # to the binary. A curve without its commit cannot be attributed later.
  [ -f "${ICE_PERF_BIN:-$ROOT/bin/ICE-mpi}.commit" ] \
    && echo "COMMIT $(tr '\n' ' ' < "${ICE_PERF_BIN:-$ROOT/bin/ICE-mpi}.commit")"
  echo "CASE  $(ls -d $CASE/INPUT* | tr '\n' ' ')"
  echo "INI   $INI md5=$(md5sum $INI | cut -c1-12)"
}

# ice_domains <core-list> <ranks> <threads> -- an explicit Intel MPI domain per rank.
#
# A `numactl --physcpubind` wrapped around mpirun does NOT pin ranks under Slurm:
# Hydra launches them through Slurm's bootstrap, which builds fresh tasks from
# the allocation and never inherits mpirun's affinity. Intel MPI then pins them
# itself over the whole node and the mask is silently ignored. Hand the
# placement to Intel MPI directly, as one hex bitmask per rank; the OpenMP
# threads then bind inside their rank's domain.
ice_domains () {
  python3 -c '
import sys
cores = [int(x) for x in sys.argv[1].split(",")]
R, T = int(sys.argv[2]), int(sys.argv[3])
if len(cores) != R * T:
    sys.exit("ice_domains: %d cores for %d ranks x %d threads" % (len(cores), R, T))
out = []
for r in range(R):
    m = 0
    for c in cores[r*T:(r+1)*T]:
        m |= 1 << c
    out.append("0x%x" % m)
print("[" + ",".join(out) + "]")' "$1" "$2" "$3"
}

# No stale ranks from the previous point: a run that starts while the previous
# one is still on the cores measures contention, not scaling.
ice_wait_clear () {
  local n
  for i in $(seq 1 60); do
    n=$(pgrep '[I]CE-mpi' 2>/dev/null | wc -l)
    [ "$n" -eq 0 ] && return 0
    [ "$i" -eq 1 ] && echo "IWARN stale ICE ranks=$n -- waiting"
    sleep 5
  done
  echo "IWARN stale ICE ranks=$n -- PROCEEDING ANYWAY"
}

# Where did this point ACTUALLY run? Emitted per row, because a campaign that
# does not check this can ship a whole curve measured on the wrong cores.
ice_occupancy () {
  ps -eLo psr,comm --no-headers 2>/dev/null | awk -v tag="$1" '
    $2=="ICE-mpi" { sol[$1]++; ns++ }
    END {
      printf "IOCC tag=%s solver_threads=%d solver_cores=%d\n", tag, ns, length(sol)
    }'
}

# ice_run <tag> <input-set> <ranks> <threads> <ppn> [core-mask] [reps]
# Emits one IREP line per repetition and one IOK line with the best.
#
# Every field of the IOK line comes from the repetition with the smallest
# wall/iter. The 2026-09 campaign took the minimum wall/iter but carried the
# cycles, phases and clock of the LAST repetition, so its "GHz" column paired
# numbers from different runs (AUDIT.md, "the flaw both harnesses share").
ice_run () {
  local tag=$1 set=$2 R=$3 T=$4 ppn=$5
  local mc reps
  # An empty-but-supplied mask must NOT degrade to "no mask": a mask that failed
  # to build would otherwise produce an unpinned run that still looks like a
  # result. Distinguish "not given" from "given and empty".
  if [ $# -ge 6 ]; then
    mc=$6; [ -z "$mc" ] && { echo "IFAIL tag=$tag reason=empty-core-mask"; return 1; }
  else mc=NONE; fi
  reps=${7:-1}

  local best="" best_fields=""
  for rep in $(seq 1 $reps); do
    local D=$W/work/${tag}-r${rep}
    rm -rf $D && mkdir -p $D/OUTPUT
    ln -sfn $CASE/$set $D/INPUT
    cp $INI $D/input.ini
    cd $D

    # Per-rank memory. Every rank holds the whole domain, so this is the number
    # that decides whether a rank count is possible at all. /usr/bin/time -v
    # around mpirun reports mpirun's own RSS and is useless here, so sample
    # /proc instead: the largest single rank, and the node total.
    ( maxrss=0; maxsum=0
      while :; do
        rss=$(ps -o rss= -C ICE-mpi 2>/dev/null)
        [ -z "$rss" ] || {
          hi=$(echo "$rss" | sort -n | tail -1)
          sum=$(echo "$rss" | awk '{s+=$1} END{print s+0}')
          [ "$hi"  -gt "$maxrss" ] 2>/dev/null && maxrss=$hi
          [ "$sum" -gt "$maxsum" ] 2>/dev/null && maxsum=$sum
        }
        echo "$maxrss $maxsum" > rss.txt
        sleep 5
      done ) &
    local rsspid=$!

    export OMP_NUM_THREADS=$T OMP_PLACES=cores
    if [ "$mc" = "NONE" ]; then
      export I_MPI_PIN_DOMAIN=omp
    else
      export I_MPI_PIN_DOMAIN=$(ice_domains "$mc" "$R" "$T") || {
        echo "IFAIL tag=$tag rep=$rep reason=domain-build"; kill $rsspid 2>/dev/null; continue; }
    fi
    echo "IPIN tag=$tag I_MPI_PIN_DOMAIN=$I_MPI_PIN_DOMAIN"
    ice_wait_clear

    local t0=$SECONDS
    # ICE_WRAP wraps the BINARY, not mpirun: a numactl around mpirun would set
    # the policy of the launcher, which Hydra does not pass to the ranks it
    # spawns through Slurm's bootstrap. Empty by default, so every existing
    # campaign launches exactly as before.
    mpirun -n $R -ppn $ppn ${ICE_WRAP:-} $BIN > run.log 2>&1 < /dev/null &
    local mpid=$!
    ( sleep 120; ice_occupancy "$tag" ) &
    local opid=$!
    wait $mpid; local rc=$?
    kill $opid 2>/dev/null; wait $opid 2>/dev/null
    local dt=$(( SECONDS - t0 ))
    kill $rsspid 2>/dev/null; wait $rsspid 2>/dev/null
    local rss_rank_gb=$(awk '{printf "%.2f", $1/1048576}' rss.txt 2>/dev/null)
    local rss_node_gb=$(awk '{printf "%.2f", $2/1048576}' rss.txt 2>/dev/null)

    # Keep the log: a campaign that deletes the evidence cannot be audited, and
    # 10 KB per point is nothing next to the case. fetch-results.sh brings
    # these back beside the .out file.
    mkdir -p $W/logs && cp run.log $W/logs/${tag}-r${rep}.log 2>/dev/null

    # A run that dies AFTER its last timer window still leaves a parsable
    # wall/iter, so the exit status has to be checked or a failed run gets
    # scored. The end-of-run summary is the second witness: it is the last
    # thing ICE prints, so its absence means the run did not finish.
    local done=$(grep -c 'Cells per rank-second' run.log)
    if [ "$rc" -ne 0 ] || [ "$done" -eq 0 ]; then
      echo "IFAIL tag=$tag rep=$rep rc=$rc finished=$done total_s=$dt -- log kept"
      tail -12 run.log
      cd $W; rm -rf $D
      continue
    fi

    # The LAST timer window, not the first: the first carries set-up transients.
    local last=$(grep 'ICE Timing |' run.log | tail -1)
    if [ -n "$last" ]; then
      local wps=$(echo "$last" | grep -o 'wall/iter *[0-9.Ee+-]*'       | awk '{print $2}')
      local tmn=$(echo "$last" | grep -o 'rank min *[0-9.Ee+-]*'        | awk '{print $3}')
      local tav=$(echo "$last" | grep -o 'avg *[0-9.Ee+-]*'             | awk '{print $2}')
      local imb=$(echo "$last" | grep -o 'imbalance *[0-9.Ee+-]*'       | awk '{print $2}')
      local cwt=$(echo "$last" | grep -o 'exchange wait *[0-9.Ee+-]*'   | awk '{print $3}')
      local swt=$(echo "$last" | grep -o 'collective wait *[0-9.Ee+-]*' | awk '{print $3}')

      local ph=$(grep 'ICE Phases |' run.log | tail -1)
      local src=$(echo "$ph" | sed 's/.*| source *//' | awk '{print $1}')
      local srp=$(echo "$ph" | sed 's/.*| source *//' | awk '{print $4}')
      local flx=$(echo "$ph" | sed 's/.*| flux *//'   | awk '{print $1}')
      local flp=$(echo "$ph" | sed 's/.*| flux *//'   | awk '{print $4}')
      local hal=$(echo "$ph" | sed 's/.*| halo *//'   | awk '{print $1}')
      local hap=$(echo "$ph" | sed 's/.*| halo *//'   | awk '{print $4}')

      local rk=$(grep 'ICE Ranks  |' run.log | tail -1)
      local wmx=$(echo "$rk" | grep -o 'compute/iter max *[0-9.Ee+-]*' | awk '{print $3}')
      local wmn=$(echo "$rk" | grep -o '| min *[0-9.Ee+-]*'            | awk '{print $3}')
      local wav=$(echo "$rk" | grep -o 'mean *[0-9.Ee+-]*'             | awk '{print $2}')
      local wsp=$(echo "$rk" | grep -o 'spread *[0-9.Ee+-]*'           | awk '{print $2}')

      local cyc=$(grep 'ICE Cycles' run.log | tail -1 \
                  | grep -o 'core-cycles/iter *[0-9.Ee+-]*' | awk '{print $2}')
      local ghz=$(grep 'ICE Cycles' run.log | tail -1 \
                  | grep -o 'aggregate *[0-9.]*' | awk '{print $2}')
      local bal=$(grep -m1 'MPI partition' run.log | grep -o 'balance *[0-9.]*%' | awk '{print $2}')
      local nbl=$(grep -m1 'MPI partition' run.log | awk '{print $3}')

      # The per-region detail line of the extended timers (absent on binaries
      # that predate it): carried verbatim, collect.py splits it.
      local det=$(grep 'ICE Detail |' run.log | tail -1 | sed 's/.*ICE Detail | *//' | tr -s ' ' | tr ' ' ';')

      # Set-up is everything outside the iteration loop: ICE reads the whole grid
      # on every rank, so it grows with rank count and is worth keeping. Take the
      # loop's own elapsed time from ICE rather than iterations x wall/iter -- the
      # cost per iteration is not flat (the reconstruction retries where the
      # limiter would make a negative density, and that grows as gradients
      # develop), so the product would be wrong and can even exceed the run.
      local loop=$(grep 'Loop elapsed, with I/O' run.log | tail -1 | awk '{print $(NF-1)}')
      local nit=$(grep -E '^ *Iterations' run.log | tail -1 | awk '{print $NF}')
      local setup=$(python3 -c "
try: print('%.1f' % (float('$dt') - float('$loop')))
except Exception: print('')" 2>/dev/null)

      # Every timer window, not only the last one. wall_iter is the last window
      # and stays so (every table of this campaign reads it), but a run that
      # starts slow and settles -- the 2026-10-01 start-up crawl: 2.1 s per
      # iteration for 20-30 iterations, then 0.35 -- hands the last window alone
      # a number that describes nothing. windows= lists them all; steady_iter is
      # the mean of the second half of the windows; transient_pct says by how
      # much the slowest window of the FIRST half exceeds that steady value (a
      # flat run reads ~20 %, the set-up effect of the first window; the crawl
      # read 500-1000 %). A point whose transient_pct is large is not a steady
      # state and must not enter a scaling table as one.
      local win=$(grep 'ICE Timing |' run.log | grep -o 'wall/iter *[0-9.Ee+-]*' | awk '{printf "%s%s", (NR>1?",":""), $2}')
      local steady=$(python3 -c "
w=[float(x) for x in '$win'.split(',') if x]
h=max(1, len(w)//2); s=sum(w[h:])/len(w[h:]) if len(w)>1 else w[0]
print('%.4E %.0f' % (s, 100.0*(max(w[:h])/s-1.0)))" 2>/dev/null)
      local std_iter=${steady%% *} trn_pct=${steady##* }

      echo "IREP tag=$tag rep=$rep rc=$rc wall_iter=$wps steady_iter=$std_iter transient_pct=$trn_pct iters=$nit total_s=$dt loop_s=$loop setup_s=$setup rss_rank_gb=$rss_rank_gb rss_node_gb=$rss_node_gb cycles_iter=$cyc ghz=$ghz windows=$win"
      local fields="steady_iter=$std_iter transient_pct=$trn_pct windows=$win rank_min_s=$tmn rank_avg_s=$tav imbalance_pct=$imb"
      fields="$fields exchange_wait_pct=$cwt collective_wait_pct=$swt"
      fields="$fields source_s=$src source_pct=$srp flux_s=$flx flux_pct=$flp halo_s=$hal halo_pct=$hap"
      fields="$fields compute_max_s=$wmx compute_min_s=$wmn compute_mean_s=$wav compute_spread_pct=$wsp"
      fields="$fields cycles_iter=$cyc ghz=$ghz blocks=$nbl balance=$bal"
      fields="$fields rss_rank_gb=$rss_rank_gb rss_node_gb=$rss_node_gb setup_s=$setup loop_s=$loop iters=$nit rep=$rep"
      [ -n "$det" ] && fields="$fields detail=$det"
      if [ -z "$best" ] || python3 -c "import sys; sys.exit(0 if float('$wps') < float('$best') else 1)"; then
        best=$wps; best_fields=$fields
      fi
    else
      echo "IFAIL tag=$tag rep=$rep total_s=$dt"; tail -8 run.log
    fi
    cd $W; rm -rf $D
  done
  if [ -n "$best" ]; then
    echo "IOK tag=$tag set=$set ranks=$R threads=$T ppn=$ppn reps=$reps wall_iter=$best $best_fields"
  fi
}
