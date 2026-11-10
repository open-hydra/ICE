#!/bin/bash

set -e  # Exit on any command failure
set -u  # Treat unset variables as an error

PROGRAM=$(basename "$0")
readonly DIR=$(pwd)
VERBOSE=false
BUILD_DIR="$DIR/build"
project=ICE

function usage() {
    cat <<EOF

Install script for $project

Usage:
  $PROGRAM [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  -h       , --help         Show this help message and exit
  -v       , --verbose      Enable verbose output

Commands:
  build                     Perform a full build
    --include-orion=<path>  Build an external ORION clone instead of the release ICE selects
    --include-finer=<path>  Build an external FiNeR clone instead of the submodule lib/FiNeR
    --update-deps           Select again the newest release of ORION (deps.lock)
    --compilers=<name>      Set compilers (intel, gnu)
    --use-openmp            Use OpenMP
    --use-mpi               Use MPI
    --use-tecio             Use TecIO

  compile                   Compile the program using the CMakePresets file

EOF
    exit 1
}

log() {
    if [ "$VERBOSE" = true ]; then
        # Bold and dim gray (ANSI escape: bold + color 90)
        echo -e "\033[1;90m$1\033[0m"
    fi
}

error() {
    # Bold red + [ERROR] tag, output to stderr
    echo -e "\033[1;31m[ERROR] $1\033[0m" >&2
}

task() {
    # Bold yellow + ==> tag, output to stdout
    echo -e "\033[1;38;5;186m==> $1\033[0m"
}


# Create default CMakePresets.json if it doesn't exist
function write_presets() {
  # the type is FILEPATH when CMake finds the compiler, STRING when it is given with -D (check_language: another)
  FC=$(grep -E '^CMAKE_Fortran_COMPILER:[A-Z]+=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CC=$(grep -E '^CMAKE_C_COMPILER:[A-Z]+=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CXX=$(grep -E '^CMAKE_CXX_COMPILER:[A-Z]+=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  # the external clones, when given
  local sources=""
  if [[ -n $ORION_PATH ]]; then
    sources+="        \"FETCHCONTENT_SOURCE_DIR_ORION\": \"${ORION_PATH}\","$'\n'
  fi
  if [[ -n $FINER_PATH ]]; then
    sources+="        \"FETCHCONTENT_SOURCE_DIR_FINER\": \"${FINER_PATH}\","$'\n'
  fi

  cat <<EOF > CMakePresets.json
{
  "version": 3,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 23
  },
  "configurePresets": [
    {
      "name": "default",
      "description": "Default preset",
      "binaryDir": "\${sourceDir}/build",
      "cacheVariables": {
${sources}        "CMAKE_BUILD_TYPE": "${BUILD_TYPE}",
        "CMAKE_Fortran_COMPILER": "${FC}",
        "CMAKE_C_COMPILER": "${CC}",
        "CMAKE_CXX_COMPILER": "${CXX}",
        "USE_TECIO": "${USE_TECIO}",
        "USE_OPENMP": "${USE_OPENMP}",
        "USE_MPI": "${USE_MPI}"
      }
    }
  ]
}
EOF
}


# Default global values
COMMAND=""
COMPILERS=""
# ORION and FiNeR: the release ICE selects and the submodule (see CMakeLists.txt), or external clones
ORION_PATH=""
FINER_PATH=""
UPDATE_DEPS="false"
USE_OPENMP="false"
USE_MPI="false"
USE_TECIO="false"
BUILD_TYPE="RELEASE"

# Define allowed options for each command using regular arrays
CMD=("build" "compile")
CMD_OPTIONS_build=("--compilers --include-orion --include-finer --update-deps --use-openmp --use-mpi --use-tecio")

# Parse global options
while getopts "hv-:" opt; do
    case "$opt" in
        -)
            case "$OPTARG" in
                verbose) VERBOSE=true ;;
                help) usage ;;
                *) error "Unknown global option '--$OPTARG'"; usage ;;
            esac
            ;;
        h) usage ;;
        v) VERBOSE=true ;;
        ?) error "Unknown global option '-$OPTARG'"; usage ;;
    esac
done
shift $((OPTIND -1))

# Ensure a command was provided
if [[ $# -eq 0 ]]; then
    error "No command provided!"
    usage
fi

COMMAND="$1"
# Check if the command is valid
if [[ ! " ${CMD[@]} " =~ " ${COMMAND} " ]]; then
    error "Unknown command '$COMMAND'"
    usage
fi
shift

# Parse command-specific options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --include-orion=*)
            [[ "$COMMAND" == "build" ]] || { error " --include-orion is only valid for 'build' command"; exit 1; }
            # absolute: CMake takes the folder from the build directory
            ORION_PATH=$(cd "${1#*=}" && pwd) || { error "No folder ${1#*=}"; exit 1; }
            ;;
        --include-finer=*)
            [[ "$COMMAND" == "build" ]] || { error " --include-finer is only valid for 'build' command"; exit 1; }
            # absolute: CMake takes the folder from the build directory
            FINER_PATH=$(cd "${1#*=}" && pwd) || { error "No folder ${1#*=}"; exit 1; }
            ;;
        --update-deps)
            [[ "$COMMAND" == "build" ]] || { error " --update-deps is only valid for 'build' command"; exit 1; }
            UPDATE_DEPS="true"
            ;;
        --compilers=*)
            [[ "$COMMAND" == "build" ]] || { error " --compilers is only valid for 'build' command"; exit 1; }
            if [[ ! "$1" =~ ^--compilers=(intel|gnu)$ ]]; then
                error "Invalid value for --compilers. Valid values are 'intel' or 'gnu'."
                exit 1
            fi
            COMPILERS="${1#*=}"
            ;;
        --use-openmp)
            [[ "$COMMAND" == "build" ]] || { error " --use-openmp is only valid for 'build' command"; exit 1; }
            USE_OPENMP="true"
            ;;
        --use-mpi)
            [[ "$COMMAND" == "build" ]] || { error " --use-mpi is only valid for 'build' command"; exit 1; }
            USE_MPI="true"
            ;;
        --use-tecio)
            [[ "$COMMAND" == "build" ]] || { error " --use-tecio is only valid for 'build' command"; exit 1; }
            USE_TECIO="true"
            ;;
        *)
            eval "opts=(\"\${CMD_OPTIONS_${COMMAND}[@]}\")"
            error "Unknown option '$1' for command '$COMMAND'. Valid options: ${opts[@]}"
            usage
            exit 1
            ;;
    esac
    shift
done


# Execute the selected command
case "$COMMAND" in
    build)
        task "Building $project"

        task "Cloning submodules"
        # sync first: a clone made before FiNeR moved (lib/third_party/FiNeR to lib/FiNeR) lacks its new entry
        git submodule sync --quiet
        [[ -z $FINER_PATH ]] && git submodule update --init lib/FiNeR
        # ORION (a release, in lib/ORION) and FiNeR's own dependencies are cloned by CMake at configure time.

        task "Configuring and building $project"
        if [[ $COMPILERS == "intel" ]]; then
          export FC="ifx"
          export CC="icx"
          export CXX="icpx"
        elif [[ $COMPILERS == "gnu" ]]; then
          export FC="gfortran"
          export CC="gcc"
          export CXX="g++"
        fi
        log "Build dir: $BUILD_DIR"
        log "Build type: $BUILD_TYPE"
        log "ORION: ${ORION_PATH:-release}"
        log "FiNeR: ${FINER_PATH:-lib/FiNeR}"
        log "Update ORION: $UPDATE_DEPS"
        log "Use OpenMP: $USE_OPENMP"
        log "Use MPI: $USE_MPI"
        log "Use TecIO: $USE_TECIO"
        if [[ -z "${FC+x}" || -z "${CXX+x}" || -z "${CC+x}" ]]; then
          log "Compilers not set. CMake will decide."
        else
          log "Compilers: FC=$FC, CXX=$CXX, CC=$CC"
        fi
        rm -rf $BUILD_DIR
        CMAKE_ARGS=(-B "$BUILD_DIR" -DUSE_TECIO="$USE_TECIO" -DUSE_OPENMP="$USE_OPENMP" -DUSE_MPI="$USE_MPI" -DCMAKE_BUILD_TYPE="$BUILD_TYPE")
        if [[ -n $ORION_PATH ]]; then
          CMAKE_ARGS+=(-DFETCHCONTENT_SOURCE_DIR_ORION="$ORION_PATH")
        fi
        if [[ -n $FINER_PATH ]]; then
          CMAKE_ARGS+=(-DFETCHCONTENT_SOURCE_DIR_FINER="$FINER_PATH")
        fi
        if [[ $UPDATE_DEPS == true ]]; then
          CMAKE_ARGS+=(-DSEMVER_UPDATE=ALL)
        fi
        cmake "${CMAKE_ARGS[@]}" || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"

        task "Write CMakePresets.json"
        write_presets
        log "[OK] CMakePresets.json created"
        ;;
    compile)
        task "Compiling $project using CMakePresets"
        cmake --preset default || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"
        ;;
    *)
        error "Unknown command '$COMMAND'"
        usage
        ;;
esac
