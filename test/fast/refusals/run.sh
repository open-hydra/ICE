#!/usr/bin/env bash
#===============================================================================
#  Fast test: a run that cannot proceed must stop with a NON-ZERO status.
#
#  Every case below is a copy of a shipped case with one thing broken, and each
#  must terminate with exit status != 0 and print the texts listed for it. The
#  first rows are refused while the input is read (a value outside its allowed
#  list, a required model left unset in a coupled run); the last one diverges at
#  run time (a NaN source in one cell) and must be caught after the update. A
#  solver that reports success on any of these would let a harness read a
#  broken run as a pass.
#===============================================================================
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../../.." && pwd)
ICE_BIN=${ICE_BIN:-$ROOT/bin/ICE}
export KMP_STACKSIZE=${KMP_STACKSIZE:-100M}

if [[ ! -x $ICE_BIN ]]; then
  echo "[fast] no ICE binary at $ICE_BIN -- build first"; exit 1
fi

WORK=$HERE/work
rm -rf "$WORK"; mkdir -p "$WORK"
fail=0

# refuse <case> <label> <mutation> <text>[;<text>...]
refuse() {
  local case=$1 label=$2 mutation=$3 texts=$4 dir rc ok t
  dir=$WORK/$(echo "$label" | tr -c 'A-Za-z0-9\n' '_')
  rm -rf "$dir"; mkdir -p "$dir"
  cp -rL "$ROOT/test/$case/INPUT" "$dir/"
  cp "$ROOT/test/$case/input.ini" "$dir/"
  ( cd "$dir" && eval "$mutation" ) || { echo "[fast] FAIL: $label (could not apply the mutation)"; fail=1; return; }
  ( cd "$dir" && mkdir -p OUTPUT && ulimit -s unlimited && OMP_NUM_THREADS=${ICE_FAST_THREADS:-1} \
    timeout 600 "$ICE_BIN" > log 2>&1 ); rc=$?
  ok=1
  [[ $rc -ne 0 ]] || ok=0
  IFS=';' read -ra need <<< "$texts"
  for t in "${need[@]}"; do grep -qF -- "$t" "$dir/log" || ok=0; done
  if [[ $ok == 1 ]]; then
    echo "[fast] PASS: $label (exit $rc)"
  else
    echo "[fast] FAIL: $label (exit $rc); expected a non-zero exit and: $texts"
    tail -5 "$dir/log" | sed 's/^/        /'
    fail=1
  fi
}

refuse Doisneau/MK "unknown flux limiter"        "sed -i 's/^flux-limiter .*/flux-limiter = Bogus/' input.ini" \
       "flux-limiter must be one of;superbee"
refuse Refuse/MK   "unknown drag model"          "sed -i 's/^drag .*/drag = Bogus/' input.ini" \
       "drag must be one of;Schiller-Naumann"
refuse Refuse/MK   "drag not set (coupled run)"  "sed -i '/^drag /d' input.ini" \
       "drag is not set;Schiller-Naumann"
refuse Refuse/MK   "heat not set (coupled run)"  "sed -i '/^heat-transfer /d' input.ini" \
       "heat-transfer is not set;Kavanau-Drake"
refuse Refuse/MK   "divergence caught after the update" ":" \
       "invalid state after the update"

exit $fail
