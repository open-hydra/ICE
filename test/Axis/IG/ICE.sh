#!/bin/bash -
#===============================================================================
#
#          FILE: ICE.sh
#
#         USAGE: run "./ICE.sh [options]" in the current shell
#
#   DESCRIPTION: A script to compile and run ICE
#===============================================================================

function print_usage {
  echo
  echo "Tasks"
  echo " compile <build>       -->     compile accordingly with the <build>"
  echo " solve                 -->     run ICE"
  echo " kill                  -->     kill the process"
  echo
  echo "Solver options"
  echo " -b | --background     -->     launch solver in background"
  echo " -m | --mpi <n>        -->     launch solver with <n> MPI ranks"
  echo " -p | --parallel <n>   -->     launch solver with <n> OMP threads per MPI rank"
  echo
  exit 1
}

# Directories and files definition
MASTERDIR=../../../
MASTER=$MASTERDIR/bin/ICE
LOCAL=./bin/ICE

# Default Options
BG=0
NMPI=1
NTHREADS=1

# Parse command-line options
while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    shift
    break
  fi
  case $1 in

    -b | --background)
        BG=1
        shift
        ;;
    -m | --mpi)
        NMPI=$2
        shift 2
        ;;
    -p | --parallel)
        NTHREADS=$2
        shift 2
        ;;
    -h | --help)
        print_usage
        shift
        ;;
    * )
      break
      ;;
  esac
done
[[ $# == 0 ]] && print_usage

DIR=$(pwd)

if [[ $1 == compile ]]; then
  mkdir -p bin
  rm -f $LOCAL
  cd $MASTERDIR
  ./install.sh compile
  cd $DIR
  cp $MASTER $LOCAL
fi

if [[ $1 == solve ]]; then
  mkdir -p OUTPUT bin
  ulimit -s unlimited
  export KMP_STACKSIZE=100M
  export OMP_NUM_THREADS=$NTHREADS
  # Check executable
  if [[ "$MASTER" -nt "$LOCAL" ]]; then
    cp $MASTER $LOCAL
  fi
  # Run the solver
  if [[ $NMPI == 1 ]]; then
    if [[ $BG == 0 ]]; then
      # keep the log: verify.py reads the boundary census and the end-of-run line from it
      $LOCAL | tee logfile
    else
      $LOCAL 2>errors_file >logfile &
      echo $! > .ID
    fi
  else
    if [[ $BG == 0 ]]; then
      mpirun -np $NMPI --map-by socket --bind-to socket $LOCAL | tee logfile
    else
      mpirun -np $NMPI --map-by socket --bind-to socket $LOCAL 2>errors_file >logfile &
      echo $! > .ID
    fi
  fi
fi

if [[ $1 == kill ]]; then
  read PID < .ID && kill $PID
fi
