#!/usr/bin/env bash
#===============================================================================
#  Bounds-checked leg: short runs on a build that checks every array subscript.
#
#  DEBUG build only (registered as Bounds, label debug, in a DEBUG configure;
#  not in CI). Five iterations each of shipped cases (one block, four blocks
#  joined by connections, two chimera blocks, an axisymmetric AG wedge, two
#  threads, and two thread groups) and of two generated boxes (make_box.py: a
#  solidifying mush, and a P2T drop on both sides of its ZGR boiling
#  temperature, whose probes must show both models at work). A run passes when
#  it exits 0, reaches its fifth iteration, prints no subscript or runtime error
#  (ifx: "severe (408)", gfortran: "Fortran runtime error"), and its banner
#  lists the bounds flag among the compiler options: the pattern
#  ICE_BOUNDS_FLAG, by default "check[ =]bound" (ifx prints "-check bound",
#  gfortran "-fcheck=bound"), so that a binary built without the checks fails
#  instead of passing for want of them.
#
#  The run on two thread groups goes with gfortran's flag (-fcheck) and is
#  skipped with ifx's: ifx 2023.1 stops it at the first step with "severe
#  (408): Subscript #1 of the array BLK has value 0" on the array-section
#  arguments of the call in compute_residual_grouped (Lib_Residual.f90), a
#  false positive: a guard `if (b < 1) error stop` put just before that call
#  never fires and the same report moves to the call's new line, and
#  gfortran's bounds-checked build runs the case clean. ICE_BOUNDS_GROUPED=1
#  runs it whatever the compiler, ICE_BOUNDS_GROUPED=0 skips it.
#===============================================================================
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
ICE_BIN=${ICE_BIN:-$ROOT/bin/ICE}
FLAG=${ICE_BOUNDS_FLAG:-check[ =]bound}
export KMP_STACKSIZE=${KMP_STACKSIZE:-100M} OMP_STACKSIZE=${OMP_STACKSIZE:-100M}

if [[ ! -x $ICE_BIN ]]; then
  echo "[fast] no ICE binary at $ICE_BIN -- build first"; exit 1
fi

WORK=$HERE/work
rm -rf "$WORK"; mkdir -p "$WORK"
fail=0

# prepare <dir> <case>: a copy of a shipped case, stopped after five iterations
prepare() {
  mkdir -p "$1/OUTPUT"
  cp -rL "$ROOT/test/$2/INPUT" "$1/"
  sed -e 's/^iter-threshold.*/iter-threshold = 5/' "$ROOT/test/$2/input.ini" > "$1/input.ini"
}

# leg <label> <dir> <threads> <groups> [<text>]: the run in <dir> must exit 0, reach its fifth iteration, print no
# subscript or runtime error, show the bounds flag in its banner, and print <text> when given
leg() {
  local label=$1 dir=$2 threads=$3 groups=$4 text=${5:-} rc why=()
  ( cd "$dir" && ulimit -s unlimited && OMP_NUM_THREADS=$threads ICE_THREAD_GROUPS=$groups \
    timeout 600 "$ICE_BIN" > log 2>&1 ); rc=$?
  [[ $rc -eq 0 ]] || why+=("exit $rc")
  grep -qE 'Iter = +5 \|' "$dir/log" || why+=("no line of iteration 5")
  ! grep -qE 'severe \(408\)|Fortran runtime error' "$dir/log" || why+=("a subscript or runtime error")
  grep -q -- "Compiled with -->.*$FLAG" "$dir/log" || why+=("no '$FLAG' among the compiler options of the banner")
  [[ -z $text ]] || grep -qF -- "$text" "$dir/log" || why+=("no '$text'")
  if [[ ${#why[@]} -eq 0 ]]; then
    echo "[fast] PASS: $label (exit $rc, OMP_NUM_THREADS=$threads, ICE_THREAD_GROUPS=$groups)"
  else
    echo "[fast] FAIL: $label: $(IFS=';'; echo "${why[*]}")"
    grep -E -A3 'severe \(408\)|Fortran runtime error' "$dir/log" | head -8 | sed 's/^/        /'
    tail -3 "$dir/log" | sed 's/^/        /'
    fail=1
  fi
}

for c in NoExchange/MK Doisneau/MK Doisneau/IG-split4 Doisneau/IG-chimera Axis/AG; do
  d=$WORK/$(echo "$c" | tr '/' '_')
  prepare "$d" "$c"
  leg "$c" "$d" 2 1
done

grouped=${ICE_BOUNDS_GROUPED:-}
if [[ -z $grouped ]]; then
  grouped=0
  grep -q -- 'Compiled with -->.* -fcheck' "$WORK/NoExchange_MK/log" && grouped=1
fi
if [[ $grouped == 1 ]]; then
  d=$WORK/Doisneau_IG-split4_groups2
  prepare "$d" Doisneau/IG-split4
  leg "Doisneau/IG-split4, two thread groups" "$d" 2 2 'OpenMP: thread groups 2 (asked 2) over 2 threads'
else
  echo "[fast] SKIP: Doisneau/IG-split4 on two thread groups (ifx 2023.1: a false severe (408) in"
  echo "       compute_residual_grouped; ICE_BOUNDS_GROUPED=1 runs it)"
fi

python3 -B "$HERE/make_box.py" solid "$WORK/solid" || { echo "[fast] FAIL: make_box.py solid"; exit 1; }
leg "solidifying mush (N4)" "$WORK/solid" 1 1

python3 -B "$HERE/make_box.py" drop "$WORK/drop" || { echo "[fast] FAIL: make_box.py drop"; exit 1; }
leg "P2T drop across its ZGR boiling temperature (M2)" "$WORK/drop" 1 1
if python3 -B "$HERE/make_box.py" check "$WORK/drop"; then
  echo "[fast] PASS: P2T drop: the ZGR and P2T surface temperatures in its probes"
else
  echo "[fast] FAIL: P2T drop: its probes do not show both models"; fail=1
fi

if [[ $fail -ne 0 ]]; then
  echo "[fast] Bounds: FAILED (runs under $WORK)"; exit 1
fi
echo "[fast] Bounds: every run within bounds"
