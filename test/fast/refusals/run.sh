#!/usr/bin/env bash
#===============================================================================
#  Fast test: a run that cannot proceed must stop with a NON-ZERO status.
#
#  Every case below is a copy of a shipped case with one thing broken, and each
#  must terminate with exit status != 0 and print the texts listed for it. The
#  first rows are refused while the input is read (a value outside its allowed
#  list, a required model left unset in a coupled run, the heat name Chang, whose
#  formula is JAXA3, a property table that is malformed or contradicts the INI);
#  the last one diverges at run time (a NaN source in one cell) and must be
#  caught after the update. A solver that reports success on
#  any of these would let a harness read a broken run as a pass. One row goes the
#  other way: an INI value equal to the table's must be accepted.
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

# accept <case> <label> <mutation> <text>[;<text>...]: the run must succeed and print the texts
accept() {
  local case=$1 label=$2 mutation=$3 texts=$4 dir rc ok t
  dir=$WORK/$(echo "$label" | tr -c 'A-Za-z0-9\n' '_')
  rm -rf "$dir"; mkdir -p "$dir"
  cp -rL "$ROOT/test/$case/INPUT" "$dir/"
  cp "$ROOT/test/$case/input.ini" "$dir/"
  ( cd "$dir" && eval "$mutation" ) || { echo "[fast] FAIL: $label (could not apply the mutation)"; fail=1; return; }
  ( cd "$dir" && mkdir -p OUTPUT && ulimit -s unlimited && OMP_NUM_THREADS=${ICE_FAST_THREADS:-1} \
    timeout 600 "$ICE_BIN" > log 2>&1 ); rc=$?
  ok=1
  [[ $rc -eq 0 ]] || ok=0
  IFS=';' read -ra need <<< "$texts"
  for t in "${need[@]}"; do grep -qF -- "$t" "$dir/log" || ok=0; done
  if [[ $ok == 1 ]]; then
    echo "[fast] PASS: $label (exit $rc)"
  else
    echo "[fast] FAIL: $label (exit $rc); expected exit 0 and: $texts"
    tail -5 "$dir/log" | sed 's/^/        /'
    fail=1
  fi
}

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
refuse Refuse/MK   "heat Chang points to JAXA3"  "sed -i 's/^heat-transfer .*/heat-transfer = Chang/' input.ini" \
       "Chang is now JAXA3;Mach-corrected law is JAXA4;- JAXA4"
# The property table (Doisneau/MK ships one, cp 1500, density 2000, and sets neither key in the INI)
refuse Doisneau/MK "table: two zones for one material" "sed -n '3,\$p' INPUT/part-properties.dat > t && cat t >> INPUT/part-properties.dat && rm t" \
       "one zone expected;expected: VARIABLES"
refuse Doisneau/MK "table: no VARIABLES line" "sed -i '/VARIABLES/d' INPUT/part-properties.dat" \
       "no VARIABLES line"
refuse Doisneau/MK "table: both enthalpy columns" "sed -i 's/\"Enthalpy\"/\"Enthalpy\", \"Enthalpy_abs\"/' INPUT/part-properties.dat" \
       "the enthalpy datum is ambiguous"
refuse Doisneau/MK "table: a short row" "sed -i '104s/ [^ ]*$//' INPUT/part-properties.dat" \
       "unreadable rows"
refuse Doisneau/MK "table: relative enthalpy with an offset" "awk 'NR<=4{print;next}{\$4=\$4+1e5;print}' INPUT/part-properties.dat > t && mv t INPUT/part-properties.dat" \
       "relative \"Enthalpy\" column with an offset"
refuse Doisneau/MK "table: rows in degrees Celsius" "sed -i 's/\"Enthalpy\"/\"Enthalpy_abs\"/' INPUT/part-properties.dat && awk 'NR<=4{print;next}{\$1=\$1-273;\$4=\$4-409500;print}' INPUT/part-properties.dat > t && mv t INPUT/part-properties.dat" \
       "a temperature below 0 K"
refuse Doisneau/MK "table: rows half a kelvin off the nodes" "awk 'NR<=4{print;next}{\$1=\$1+0.5;print}' INPUT/part-properties.dat > t && mv t INPUT/part-properties.dat" \
       "rows not on consecutive integer kelvins"
refuse Doisneau/MK "table: INI density against a constant column" "sed -i '/^heat-transfer/a density = 2700' input.ini" \
       "density = 2.70000E+03 differs from the table"
refuse Doisneau/MK "table: INI cp against a varying column" \
       "python3 -c \"L=open('INPUT/part-properties.dat').read().split(chr(10))[:4]; L+=['%.1f %.10g 2000.0 %.10g' % (T, 1500+0.1*T, 1500*T+0.05*T*T) for T in range(1, 5001)]; open('INPUT/part-properties.dat','w').write(chr(10).join(L)+chr(10))\" && sed -i '/^heat-transfer/a specific-heat = 1500' input.ini" \
       "specific-heat is set but the table"
refuse Doisneau/MK "table: text after the last row" "echo END >> INPUT/part-properties.dat" \
       "text after the last data row"
refuse Doisneau/MK "table: fewer rows announced than held" "sed -i 's/^I=5000,/I=4999,/' INPUT/part-properties.dat" \
       "the zone announces 4999 rows but holds 5000"
refuse Doisneau/MK "table: more rows announced than held" "sed -i 's/^I=5000,/I=5001,/' INPUT/part-properties.dat" \
       "the zone announces 5001 rows but holds 5000"
refuse Doisneau/MK "table: a value that is not a number" "sed -i '1000s/ [^ ]*\$/ abc/' INPUT/part-properties.dat" \
       "unreadable rows: line 1000 does not hold 4 numbers"
accept Doisneau/MK "table: INI density equal to the constant column" \
       "sed -i '/^heat-transfer/a density = 2000' input.ini && sed -i 's/^iter-threshold .*/iter-threshold = 20/' input.ini" \
       "material 1: density 2.00000E+03 (constant)"
refuse Refuse/MK   "divergence caught after the update" ":" \
       "invalid state after the update"

exit $fail
