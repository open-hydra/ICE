#===============================================================================
#  Helpers shared by the fast tests: build short runs of an existing case in a
#  scratch directory, run them, and compare their output byte for byte.
#
#  The cases themselves are the ones under test/Doisneau; only the iteration
#  count is shortened, so the fast tier needs no data of its own.
#===============================================================================

ICE_BIN=${ICE_BIN:-$ROOT/bin/ICE}
export OMP_STACKSIZE=${OMP_STACKSIZE:-100M}
export KMP_STACKSIZE=${KMP_STACKSIZE:-100M}

# ice_fast_prepare <workdir> <case> <iterations> <run-id>...
# ICE_FAST_SD=<detector> replaces the case's shock-detector key (ice_fast_detector checks the runs took it).
ice_fast_prepare() {
  local work=$1 case=$2 iter=$3; shift 3
  local src=$ROOT/test/$case id

  if [[ ! -x $ICE_BIN ]]; then
    echo "[fast] no ICE binary at $ICE_BIN -- build first"
    exit 1
  fi

  rm -rf "$work"
  for id in "$@"; do
    mkdir -p "$work/run$id/OUTPUT"
    cp -rL "$src/INPUT" "$work/run$id/"
    sed -e "s/^iter-threshold.*/iter-threshold = $iter/" \
        -e "s/^shell-diter.*/shell-diter = $iter/" \
        ${ICE_FAST_SD:+-e "s/^shock-detector.*/shock-detector = $ICE_FAST_SD/"} "$src/input.ini" > "$work/run$id/input.ini"
  done
}

# ice_fast_detector <dir>... -- with ICE_FAST_SD set (not none), every run's banner must report that detector
ice_fast_detector() {
  local d
  [[ -z ${ICE_FAST_SD:-} || ${ICE_FAST_SD} == none ]] && return 0
  for d in "$@"; do
    if ! grep -q "Shock *--> *$ICE_FAST_SD detector" "$d/log"; then
      echo "[fast] FAIL: ICE_FAST_SD=$ICE_FAST_SD, but the banner of $d/log reports no $ICE_FAST_SD detector"
      exit 1
    fi
  done
}

# ice_fast_run <dir> <threads> <ranks>
ice_fast_run() {
  local dir=$1 threads=$2 ranks=$3
  (
    cd "$dir"
    ulimit -s unlimited || true
    export OMP_NUM_THREADS=$threads
    if [[ $ranks -gt 1 ]]; then
      mpirun -np "$ranks" "$ICE_BIN" > log 2> err
    else
      "$ICE_BIN" > log 2> err
    fi
  )
}

# ice_fast_compare <dir-a> <dir-b> <label-a> <label-b>
ice_fast_compare() {
  local a=$1 b=$2 la=$3 lb=$4
  local f=OUTPUT/part-field.tec

  if [[ ! -s $a/$f || ! -s $b/$f ]]; then
    echo "[fast] FAIL: the solver produced no $f (see $a/log, $b/log)"
    exit 1
  fi
  if cmp -s "$a/$f" "$b/$f"; then
    echo "[fast] PASS: $la and $lb are bit-identical"
  else
    echo "[fast] FAIL: $la and $lb differ"
    cmp "$a/$f" "$b/$f" | head -5
    exit 1
  fi
}
