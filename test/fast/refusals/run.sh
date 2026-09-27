#!/usr/bin/env bash
#===============================================================================
#  Fast test: a run that cannot proceed must stop with a NON-ZERO status.
#
#  Every case below is a copy of a shipped case with one thing broken, and each
#  must terminate with exit status != 0 and print the texts listed for it. The
#  first rows are refused while the input is read (a value outside its allowed
#  list, a required model left unset in a coupled run, the heat name Chang, whose
#  formula is JAXA3, a property table that is malformed or contradicts the INI,
#  a Psat column an evaporation model cannot use, a key given together with its
#  alias, a phase file whose materials or model tokens ICE cannot honour);
#  the last one diverges at run time (a NaN source in one cell) and must be
#  caught after the update. A solver that reports success on
#  any of these would let a harness read a broken run as a pass. Four rows go the
#  other way: an INI value equal to the table's must be accepted, the boiling
#  temperature must be read under either of its names, and a solidifying material
#  must be read, reported and written with its frozen and nucleated fractions.
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
# Refuse/MK and NoExchange/MK are coupled, so their evaporation model is read. psat_table writes
# their material (cp 1000, density 2000, as their INIs set) as a table with a Psat column: 1e6 - T Pa
# ("decreasing"), T Pa ("linear", 0.03 atm at the default boiling temperature, 2792 K) or
# 101.325 T Pa ("1atm", one atmosphere at 1000 K)
psat_table() {
  awk -v mode="$1" 'BEGIN {
    print "TITLE = \"Mass Thermodynamic Properties\""
    print "VARIABLES = \"Temperature\", \"Cp\", \"Density\", \"Enthalpy\", \"Psat\""
    print "ZONE T=\"A\""; print "I=5000, F=POINT"
    for (T = 1; T <= 5000; T++) printf "%.1f 1000.0 2000.0 %.1f %.6e\n", T, 1000*T, \
      (mode == "decreasing" ? 1e6 - T : (mode == "1atm" ? 101.325*T : T))
  }' > INPUT/part-properties.dat
}
evap_on="sed -i '/^heat-transfer/a evaporation = CEM' input.ini"
refuse Refuse/MK   "Psat: decreasing, with evaporation" "psat_table decreasing && $evap_on" \
       "Psat column: a pressure that decreases with T"
refuse Refuse/MK   "Psat: far from 1 atm at the boiling temperature" "psat_table linear && $evap_on" \
       "Psat column: psat(boiling-temperature) is not within a factor 2 of one atmosphere"
# Tboil is the alias of boiling-temperature (IGLOO's name): one or the other, even with equal values.
# Without the alias Tboil is an unknown key: the run stops, but not with this text.
refuse Doisneau/MK "boiling temperature: key and alias both given" \
       "sed -i '/^heat-transfer/a boiling-temperature = 373.15\nTboil = 373.15' input.ini" \
       "give boiling-temperature or its alias Tboil, not both"
# The alias is read, not only accepted: the "1atm" column gives 1.2 atm at 1200 K under either name,
# and at the default 2792 K (2.792 atm) the run would be refused.
accept NoExchange/MK "boiling temperature: key read" \
       "psat_table 1atm && $evap_on && sed -i '/^heat-transfer/a boiling-temperature = 1200' input.ini" \
       "psat(boiling-temperature) = 1.20000E+00 atm"
accept NoExchange/MK "boiling temperature: alias Tboil read" \
       "psat_table 1atm && $evap_on && sed -i '/^heat-transfer/a Tboil = 1200' input.ini" \
       "psat(boiling-temperature) = 1.20000E+00 atm"
# Materials come from the phase file: "<name> <groups> [key=value ...]" per line. two_mat declares
# a second material and a second family; the tokens are read for every material, coupled or not.
two_mat="printf 'condensed-dispersed phase\nA 1\nB 1\n' > INPUT/part-phase.txt && printf '\n[ICE-Family2]\nclosure = MK\n' >> input.ini"
tokens() { printf 'condensed-dispersed phase\nA 1 %s\n' "$1" > INPUT/part-phase.txt; }
refuse Doisneau/MK "materials: more populations than families" \
       "printf 'condensed-dispersed phase\nA 1\nB 1\n' > INPUT/part-phase.txt" \
       "declares 2 populations over 2 materials, but [ICE-Family*] defines 1 families"
refuse Doisneau/MK "materials: two materials without a table" "$two_mat && rm INPUT/part-properties.dat" \
       "absent, and 2 materials need it, one zone each"
refuse Doisneau/MK "materials: one table zone for two materials" "$two_mat" \
       "2 zones expected (one per material)"
refuse Doisneau/MK "materials: a vector of the wrong length" \
       "sed -i '/^heat-transfer/a latent-heat = 1e6 2e6' input.ini" \
       "latent-heat carries 2 values for 1 material(s)"
refuse Doisneau/MK "tokens: an unknown evaporation model" "tokens evaporation=LEB" \
       "Wrong evaporation model input ---> LEB;Choose one of the following"
refuse Doisneau/MK "tokens: combustion is not modelled" "tokens combustion=Beckstead" \
       "Wrong combustion input ---> Beckstead;- none"
# Solidification, from the tokens as IGLOO reads them; "valid" is the solid-box material's line
valid='solidification=on h-fus=1.07e6 cp-solid=600'
accept Doisneau/MK "solidification: read and reported" "tokens '$valid' && sed -i 's/^iter-threshold .*/iter-threshold = 20/' input.ini" \
       " - Solid   --> material 1 (A)"
solid_out=$WORK/$(echo "solidification: read and reported" | tr -c 'A-Za-z0-9\n' '_')/OUTPUT/part-field.tec
if grep -qF '"n_p1""f_p1""chi_p1"' "$solid_out" 2>/dev/null; then
  echo "[fast] PASS: solidification: f_p1 and chi_p1 follow n_p1 in the output"
else
  echo "[fast] FAIL: solidification: the output does not name f_p1 and chi_p1 after n_p1"; fail=1
fi
refuse Doisneau/MK "solidification: a value that is not on/off" "tokens solidification=maybe" \
       "Wrong solidification input ---> maybe;- on;- off"
refuse Doisneau/MK "solidification: heat of fusion missing" "tokens 'solidification=on cp-solid=600'" \
       "solidification=on requires h-fus > 0"
refuse Doisneau/MK "solidification: solid heat capacity missing" "tokens 'solidification=on h-fus=1.07e6'" \
       "solidification=on requires cp-solid > 0"
refuse Doisneau/MK "solidification: nucleation above melting" "tokens '$valid T-nuc=2400'" \
       "solidification=on requires T-nuc < T-melt"
refuse Doisneau/MK "solidification: with evaporation" "tokens '$valid evaporation=CEM'" \
       "solidification=on is exclusive with evaporation"
refuse Doisneau/MK "solidification: a varying cp" \
       "tokens '$valid' && python3 -c \"L=open('INPUT/part-properties.dat').read().split(chr(10))[:4]; L+=['%.1f %.10g 2000.0 %.10g' % (T, 1500+0.1*T, 1500*T+0.05*T*T) for T in range(1, 5001)]; open('INPUT/part-properties.dat','w').write(chr(10).join(L)+chr(10))\"" \
       "solidification=on requires a constant-cp material"
refuse Doisneau/MK "solidification: a varying density" \
       "tokens '$valid' && python3 -c \"L=open('INPUT/part-properties.dat').read().split(chr(10))[:4]; L+=['%.1f 1500.0 %.10g %.10g' % (T, 2000+0.1*T, 1500*T) for T in range(1, 5001)]; open('INPUT/part-properties.dat','w').write(chr(10).join(L)+chr(10))\"" \
       "solidification=on requires a constant-density material"
refuse Doisneau/MK "solidification: an IC of the wrong width" \
       "tokens '$valid' && python3 -c \"import re; L=open('INPUT/part-ic.tec').read().rstrip(chr(10)).split(chr(10)); I,J,K=[int(re.search(r'\\b%s\\s*=\\s*(\\d+)' % c, L[1]).group(1)) for c in 'IJK']; n=(I-1)*(J-1)*(K-1); open('INPUT/part-ic.tec','w').write(chr(10).join(L+L[-n:])+chr(10))\"" \
       "the initial condition holds 7 variables per cell"
refuse Doisneau/MK "tokens: only the infinite-conductivity liquid" "tokens liquid-conduction=P2T" \
       "Wrong liquid-conduction input ---> P2T;- ITC"
refuse Doisneau/MK "tokens: only the boiling clamp" "tokens boiling=ZGR" \
       "Wrong boiling input ---> ZGR;- clamp"
refuse Doisneau/MK "tokens: an unknown key" "tokens colour=blue" \
       'unknown key "colour"; the keys are'
refuse Doisneau/MK "tokens: a real key that is not a number" "tokens alpha-e=abc" \
       "alpha-e=abc is not a real number"
refuse Doisneau/MK "tokens: LK interface on the d-squared law" "tokens 'evaporation=d2-law interface=LK'" \
       "interface = LK needs a gas-side evaporation model"
refuse Refuse/MK   "divergence caught after the update" ":" \
       "invalid state after the update"

exit $fail
