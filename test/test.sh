#!/bin/bash -
#===============================================================================
#
#          FILE: test.sh
#
#         USAGE: ./test.sh [options] <command> <test>
#
#   DESCRIPTION: Run ICE validation test cases and verify results
#===============================================================================

function print_usage {
  echo "Bash script to run the ICE validation tests"
  echo
  echo "Usage:"
  echo "   ./test.sh check  all           => run and verify all test cases"
  echo "   ./test.sh check  Doisneau/MK   => run and verify a specific test"
  echo "   ./test.sh update all           => run and update all reference solutions"
  echo "   ./test.sh update Doisneau/MK   => run and update a specific reference"
  echo "   ./test.sh clean                => clean all test directories"
  echo
  echo "Options:"
  echo "   -p | --parallel <n>            => run solver with <n> cores (default: 1)"
  echo
  exit 1
}

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

TESTROOT=$(pwd)
FILE=$TESTROOT/logfile
NTHREADS=1

# Parse options
while test $# -gt 0; do
  case "$1" in
    -p | --parallel)
      NTHREADS=$2
      shift 2
      ;;
    * )
      break
      ;;
  esac
done

ALL_TESTS=(
  Doisneau/MK
  Doisneau/IG
  Doisneau/AG
)

function clean {
  for TEST in "${ALL_TESTS[@]}"; do
    cd $TESTROOT/$TEST
    rm -rf bin OUTPUT logfile errors_file .ID
    cd $TESTROOT
  done
  rm -f $FILE
}

function run_solver {
  local TEST=$1
  cd $TESTROOT/$TEST
  mkdir -p bin OUTPUT
  bash ICE.sh -p $NTHREADS solve > logfile 2>errors_file
}

function check_test {
  local TEST=$1
  run_solver $TEST
  cd $TESTROOT/$TEST
  python3 -B verify.py
  local esito=$?
  if [[ $esito == "0" ]]; then
    echo "PASS - $TEST" >> $FILE
  else
    echo "FAIL - $TEST" >> $FILE
  fi
  cd $TESTROOT
}

function update_test {
  local TEST=$1
  run_solver $TEST
  cp $TESTROOT/$TEST/OUTPUT/part-field.tec $TESTROOT/$TEST/reference/part-field.tec
  echo -e "$TEST  -->  ${GREEN}Updated${NC}"
  echo "Updated - $TEST" >> $FILE
  cd $TESTROOT
}

[[ $# == 0 ]] && print_usage

CMD=$1
TARGET=$2

if [[ $CMD == clean ]]; then
  clean
  exit 0
fi

[[ -z $TARGET ]] && print_usage

rm -f $FILE
echo 'ICE TEST SESSION' >> $FILE
date >> $FILE
echo >> $FILE

if [[ $TARGET == all ]]; then
  for TEST in "${ALL_TESTS[@]}"; do
    [[ $CMD == check  ]] && check_test  $TEST
    [[ $CMD == update ]] && update_test $TEST
  done
else
  [[ $CMD == check  ]] && check_test  $TARGET
  [[ $CMD == update ]] && update_test $TARGET
fi
