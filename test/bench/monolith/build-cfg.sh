#!/usr/bin/env bash
# Build one ICE configuration into its own build and bin directories.
#   build-cfg.sh <cfg> [<ICE root>]      cfg: serial | omp | mpi | hyb | omp-strict |
#                                             hyb-strict | omp-ipo | omp-tmp
# Works on sprop2 and on monolith; sync-src-monolith.sh calls it over ssh.
#
# Notes that cost time if rediscovered (Marco Grossi, 2026-09, and 2026-09-30):
#  * On monolith Intel oneAPI 2023.0.0 + Intel MPI are in the LOGIN environment
#    only: run through `bash -l`. On sprop2 anaconda's mpirun/mpifort shadow
#    Intel's: name the wrapper by absolute path (MPIFC) and export I_MPI_F90=ifx,
#    or FindMPI pairs mpiifort's ifort with an ifx build.
#  * TecIO is C++, and an ifx link does not pull in libstdc++ by itself. It has
#    to arrive AFTER libteciompi.a on the link line -- a static archive only
#    pulls what is already unresolved -- so it goes in
#    CMAKE_Fortran_STANDARD_LIBRARIES, which CMake appends last, not in
#    CMAKE_EXE_LINKER_FLAGS, which it puts first and which therefore silently
#    does nothing here.
#  * The C compiler is needed only for src/lib/diagnostic/ice_cycles.c (the
#    per-thread cycle counter); it is plain POSIX C and links through the C ABI.
#  * CXX must be named too. ORION installs TecIO into a directory suffixed with
#    the C++ and Fortran compiler IDs; pinning icpx keeps one install directory
#    rather than two half-populated ones. ORION builds TecIO from
#    lib/ORION/lib/TecIO/teciompisrc on first configure.
#  * ICE_BIN_DIR puts each configuration's binary in its own directory; before
#    it every configuration overwrote bin/ICE.
set -eu
CFG=${1:?configuration}
SRC=${2:-$(cd "$(dirname "$0")/../../.." && pwd)}
FC=${FC:-$(command -v ifx || command -v ifort)}
MPIFC=${MPIFC:-$(command -v mpiifort)}
CC=${CC:-$(command -v icx || command -v gcc)}
CXX=${CXX:-$(command -v icpx || command -v g++)}
export I_MPI_F90=$FC

TYPE=RELEASE; FLAGS=""
case $CFG in
  serial)     OPTS="-DUSE_MPI=OFF -DUSE_OPENMP=OFF" ;;
  omp)        OPTS="-DUSE_MPI=OFF -DUSE_OPENMP=ON" ;;
  mpi)        OPTS="-DUSE_MPI=ON  -DUSE_OPENMP=OFF" ;;
  hyb)        OPTS="-DUSE_MPI=ON  -DUSE_OPENMP=ON" ;;
  omp-strict) OPTS="-DUSE_MPI=OFF -DUSE_OPENMP=ON"; FLAGS="-fp-model=strict" ;;
  hyb-strict) OPTS="-DUSE_MPI=ON  -DUSE_OPENMP=ON"; FLAGS="-fp-model=strict" ;;
  omp-ipo)    OPTS="-DUSE_MPI=OFF -DUSE_OPENMP=ON -DICE_ENABLE_IPO=ON" ;;
  hyb-ipo)    OPTS="-DUSE_MPI=ON  -DUSE_OPENMP=ON -DICE_ENABLE_IPO=ON" ;;
  omp-tmp)    OPTS="-DUSE_MPI=OFF -DUSE_OPENMP=ON"; TYPE=DEBUG; FLAGS="-check arg_temp_created" ;;
  *) echo "unknown configuration $CFG" >&2; exit 2 ;;
esac

cmake -S "$SRC" -B "$SRC/build-$CFG" \
  -DCMAKE_BUILD_TYPE=$TYPE \
  -DCMAKE_Fortran_COMPILER="$FC" \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DMPI_Fortran_COMPILER="$MPIFC" \
  -DCMAKE_Fortran_STANDARD_LIBRARIES=-lstdc++ \
  -DICE_BIN_DIR="$SRC/bin-$CFG" \
  ${FLAGS:+-DCMAKE_Fortran_FLAGS="$FLAGS"} \
  $OPTS

cmake --build "$SRC/build-$CFG" -j "${JOBS:-16}"

echo "== $CFG: $SRC/bin-$CFG/ICE"
md5sum "$SRC/bin-$CFG/ICE"
ldd "$SRC/bin-$CFG/ICE" | grep -E 'libmpi|libiomp5' || echo "   (no libmpi/libiomp5 linked)"
