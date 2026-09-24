# ---------------------------------------------------------------------------
#  CTest definitions for the ICE test cases
#
#  Tiers are expressed through LABELS:
#    fast        -> short runs of hard invariants; run by the pre-push hook
#                   (ctest -L fast)
#    validation  -> full case compared against its stored reference solution
#    needs-mpi   -> requires an MPI build (USE_MPI=ON); registered only then
#
#  The cases run in their own source directory and leave their output in it,
#  exactly as when they are run by hand with ./ICE.sh.
# ---------------------------------------------------------------------------

set(ICE_TEST_THREADS 1 CACHE STRING "OpenMP threads per rank used by the CTest cases")
set(ICE_TEST_RANKS   1 CACHE STRING "MPI ranks used by the CTest validation cases (needs USE_MPI)")

set(_ICE_TEST_ENV "KMP_STACKSIZE=100M")

# Register a case: run the solver through its ICE.sh, then its verify.py.
function(ice_add_case)
    cmake_parse_arguments(T "" "NAME;DIR;TIMEOUT" "LABELS" ${ARGN})
    if(NOT T_TIMEOUT)
        set(T_TIMEOUT 1800)
    endif()
    add_test(
        NAME    ${T_NAME}
        COMMAND bash -c "./ICE.sh -m ${ICE_TEST_RANKS} -p ${ICE_TEST_THREADS} solve && python3 -B verify.py"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/${T_DIR}"
    )
    set_tests_properties(${T_NAME} PROPERTIES
        LABELS      "${T_LABELS}"
        TIMEOUT     ${T_TIMEOUT}
        ENVIRONMENT "${_ICE_TEST_ENV}"
    )
endfunction()

# Register a code-verification case: a Python driver that builds its own inputs,
# runs the solver and compares against a solution known independently of ICE.
function(ice_add_verification NAME DIR)
    add_test(
        NAME    ${NAME}
        COMMAND bash -c "python3 -B run.py"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/verification/${DIR}"
    )
    set_tests_properties(${NAME} PROPERTIES
        LABELS      "${ARGN}"
        TIMEOUT     900
        ENVIRONMENT "${_ICE_TEST_ENV}"
    )
endfunction()

# Register a self-contained fast test driven by its own run.sh.
function(ice_add_fast_script NAME DIR)
    add_test(
        NAME    ${NAME}
        COMMAND bash -c "./run.sh"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test/${DIR}"
    )
    set_tests_properties(${NAME} PROPERTIES
        LABELS      "${ARGN}"
        TIMEOUT     600
        ENVIRONMENT "${_ICE_TEST_ENV}"
    )
endfunction()

# ===========================================================================
#  VALIDATION tier  (ctest -L validation) - the crossing-jets suite
# ===========================================================================

# Free-streaming crossing jets, one case per closure (docs/vv/crossing-jets.md)
ice_add_case(NAME DoisneauMK DIR Doisneau/MK LABELS "validation;MK;2D")
ice_add_case(NAME DoisneauIG DIR Doisneau/IG LABELS "validation;IG;2D")
ice_add_case(NAME DoisneauAG DIR Doisneau/AG LABELS "validation;AG;2D")

# Same case on two overlapping non-matching blocks (docs/vv/chimera.md)
ice_add_case(NAME DoisneauChimera DIR Doisneau/IG-chimera LABELS "validation;IG;chimera;2D")

# Same case on four blocks joined by connections (docs/vv/multiblock-mpi.md)
ice_add_case(NAME DoisneauSplit4 DIR Doisneau/IG-split4 LABELS "validation;IG;connection;2D")

# ===========================================================================
#  VERIFICATION tier  (ctest -L verification) - exact and reference solutions
#
#  These compare ICE against solutions that do not come from ICE, so they can
#  fail on day one; the validation cases above can only detect a change.
# ===========================================================================

# Drag and heat source terms against their closed-form relaxation
ice_add_verification(DragRelaxation    A-drag-relaxation   "fast;verification;sources")
ice_add_verification(ThermalRelaxation B-thermal-relaxation "fast;verification;sources")

# Eulerian transport against free advection, with a grid-refinement study
ice_add_verification(SineAdvection     C-advection-sine    "fast;verification;transport")

# Every drag and Nusselt correlation against an independent integration. Slower
# than the fast tier allows, since each law is a separate run.
ice_add_verification(DragMatrix        D-drag-matrix       "verification;sources")
ice_add_verification(HeatMatrix        E-heat-matrix       "verification;sources")

# Every evaporation model against the d-squared law it reduces to at zero slip, then
# the coupled mass-and-temperature system against an independent integration
# (docs/vv/verification.md).
ice_add_verification(Evaporation       F-evaporation       "verification;sources;phase-change")

# A cloud released into a frozen gas: transport and drag together, against the exact
# rigid translation the two produce (docs/vv/verification.md).
ice_add_verification(CloudTranslation  G-cloud-translation "verification;sources;transport")

# The same cloud in a gas that strains it, against the exact affine solution, and its
# slip against the small-Stokes asymptote u_p ~ u_g - tau Du_g/Dt.
ice_add_verification(LinearStrain      H-linear-strain     "verification;sources;transport")

# Two-dimensional: a cloud orbiting in a prescribed solid-body vortex, against the
# exact conformal map, swept over four decades of Stokes number. The only verification
# case that exercises the j-direction fluxes.
ice_add_verification(VortexCloud       I-vortex-cloud      "verification;sources;transport;2D")

# Berthon's Riemann problems for the anisotropic Gaussian closure, against the
# analytical wave patterns (docs/vv/berthon.md). These are full case directories
# rather than self-building drivers, so they go through ice_add_case; they are the
# only fast-tier cases that exercise AG, its Rusanov solver and the MC limiter.
ice_add_case(NAME BerthonSCS DIR Berthon/SCS LABELS "fast;verification;AG;transport;1D")
ice_add_case(NAME BerthonRCS DIR Berthon/RCS LABELS "fast;verification;AG;transport;1D")
ice_add_case(NAME BerthonRCR DIR Berthon/RCR LABELS "fast;verification;AG;transport;1D")

# The same problem through the HLLE flux, which the `riemann` key selects in place of
# the solver the closure would pick. Nothing else exercises that path.
ice_add_case(NAME BerthonHLLE DIR Berthon/SCS-hlle LABELS "verification;AG;transport;1D")

# ===========================================================================
#  FAST tier  (pre-push: ctest -L fast) - short runs, no stored references
# ===========================================================================

# Thread count must not change the answer: 1 vs 4 threads, bit-identical
ice_add_fast_script(OpenMPEquiv fast/openmp-equiv "fast;implementation")

# A run that cannot proceed must stop with a non-zero status (broken inputs, a
# required model left unset in a coupled run, a NaN after the update).
ice_add_fast_script(Refusals fast/refusals "fast;implementation;refusal")

# Rank count must not change the answer either. Needs an MPI build, so it is
# registered only when there is one.
if(USE_MPI)
    ice_add_fast_script(MPIEquivConnection fast/mpi-equiv "fast;implementation;needs-mpi")
    set_tests_properties(MPIEquivConnection PROPERTIES
        ENVIRONMENT "${_ICE_TEST_ENV};ICE_FAST_CASE=Doisneau/IG-split4")

    ice_add_fast_script(MPIEquivChimera fast/mpi-equiv "fast;implementation;chimera;needs-mpi")
    set_tests_properties(MPIEquivChimera PROPERTIES
        ENVIRONMENT "${_ICE_TEST_ENV};ICE_FAST_CASE=Doisneau/IG-chimera")
endif()
