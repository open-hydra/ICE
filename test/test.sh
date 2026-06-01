#!/bin/bash -
#===============================================================================
#
#          FILE: test.sh
#
#         USAGE: ./test.sh [options]
#
#   DESCRIPTION: Run ICE validation test cases and verify results
#===============================================================================

function print_usage {
  echo "Bash script to run the ICE validation tests"
  echo
  echo "Usage:"
  echo "   ./test.sh all              => run all test cases"
  echo "   ./test.sh Doisneau/MK      => run a specific test case"
  echo "   ./test.sh clean            => clean all test directories"
  echo
  echo "Options:"
  echo "   -p | --parallel <n>        => run solver with <n> cores (default: 4)"
  echo
  exit 1
}

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

TESTROOT=$(pwd)
FILE=$TESTROOT/testlog.txt
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

function run_test {
  local TEST=$1
  cd $TESTROOT/$TEST
  mkdir -p bin OUTPUT

  bash ICE.sh -p $NTHREADS solve > logfile 2>errors_file

  python3 -B verify.py
  local esito=$?
  if [[ $esito == "0" ]]; then
    echo "PASS - $TEST" >> $FILE
  else
    echo "FAIL - $TEST" >> $FILE
  fi

  cd $TESTROOT
}

[[ $# == 0 ]] && print_usage

if [[ $1 == clean ]]; then
  clean
  exit 0
fi

rm -f $FILE
echo 'ICE TEST SESSION' >> $FILE
date >> $FILE
echo >> $FILE

if [[ $1 == all ]]; then
  for TEST in "${ALL_TESTS[@]}"; do
    run_test $TEST
  done
else
  run_test $1
fi
