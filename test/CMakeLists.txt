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

# Particles whose relaxation time is below dt-max: tau-factor bounds the step by it in
# both stepping modes, and cells below the empty-cell density do not set it.
ice_add_verification(TauLimit          J-tau-limit         "fast;verification;sources")

# The property table read at one fixed temperature through the drag relaxation: linear
# interpolation between the integer-kelvin rows, a table from Tmin > 1, saturation below it.
ice_add_verification(PropertiesTable   K-properties-table  "fast;verification;sources")

# Several materials: the phase file maps the families onto materials, each with its own
# table zone, model tokens and inlet records, against each material's own oracle.
ice_add_verification(Materials         L-materials         "fast;verification;sources")

# Solidification: a closed cell through supercooling, recalescence, the plateau, the solid and melting at T-melt, against
# the closed forms; a steady stream against its first-order recursion, the closed form and a Rusanov model, a heated
# solid melting, a gas change that must move the nucleation front back (no lock), and the stream on two blocks.
ice_add_verification(SolidCooling      N-solid-cooling     "fast;verification;sources;phase-change")
ice_add_verification(SolidAdvection    O-solid-advection   "fast;verification;sources;transport;phase-change")
set_tests_properties(SolidAdvection PROPERTIES RESOURCE_LOCK O-solid-advection-work)

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

# A coupled slab with drag = NoDrag and heat-transfer = NoHeat: the particles must keep
# their initial state bit for bit although the gas moves and is hotter.
ice_add_case(NAME NoExchangeMK DIR NoExchange/MK LABELS "fast;verification;MK;sources;2D")

# A one-degree wedge about the x axis: the ghosts of the wedge side faces (code 200) are
# mirrored, and for IG/AG the pressure on those faces balances the radial faces, so a
# cloud at rest stays at rest.
ice_add_case(NAME AxisMK DIR Axis/MK LABELS "fast;verification;MK;2D")
ice_add_case(NAME AxisIG DIR Axis/IG LABELS "fast;verification;IG;2D")
ice_add_case(NAME AxisAG DIR Axis/AG LABELS "fast;verification;AG;2D")

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

# A pressurised blob expanding from rest in a closed box at uniform temperature: with
# the pressure work in the energy flux the thermal energy moves with the mass, so the
# temperature of the Gaussian closures stays uniform to round-off.
ice_add_case(NAME ThermalIG DIR Thermal/IG LABELS "fast;verification;IG;2D")
ice_add_case(NAME ThermalAG DIR Thermal/AG LABELS "fast;verification;AG;2D")

# A strongly anisotropic Sod tube (P11/P22 = 1e7) through explicit Euler and first order,
# against the exact gamma = 3 solution the 1-D AG system reduces to. Free of oscillations
# only when Rusanov and the time step use the directional speed sqrt(3 P_nn / rho).
ice_add_case(NAME RiemannAG DIR Riemann/AG LABELS "fast;verification;AG;transport;1D")

# A sheared AG blob on a code-200 symmetry plane: the mirrored tensor exerts no shear
# on the plane, so the x-momentum is conserved to round-off.
ice_add_case(NAME ReflectAG DIR Reflect/AG LABELS "fast;verification;AG;2D")

# ===========================================================================
#  FAST tier  (pre-push: ctest -L fast) - short runs, no stored references
# ===========================================================================

# Thread count must not change the answer: 1 vs 4 threads, bit-identical
ice_add_fast_script(OpenMPEquiv fast/openmp-equiv "fast;implementation")

# Which copy of the ATLAS boundary table belongs to which particle family: one
# copy speaks for all of them, one copy per family speaks for its own, anything
# else stops the setup. Builds its own cases, so it owns no data.
ice_add_fast_script(PopulationBinding fast/populations "fast;implementation;bc")

# Grid sequencing: the multigrid levels run, and neither the thread count nor the
# rank count changes the answer. Nothing else in the suite sets MG-levels > 1, so
# this is the only cover the restriction, the prolongation, the per-level boundary
# tables and the per-level halo schedules have. Builds its own case.
ice_add_fast_script(Multigrid fast/multigrid "fast;implementation;multigrid")

# A run that cannot proceed must stop with a non-zero status (broken inputs, a
# required model left unset in a coupled run, a NaN after the update).
ice_add_fast_script(Refusals fast/refusals "fast;implementation;refusal")

# Unit tests: programs linked against the library, each run in its own directory
function(ice_add_unit_test NAME SOURCE)
    add_executable(${NAME} ${SOURCE})
    set_target_properties(${NAME} PROPERTIES LINKER_LANGUAGE Fortran
                          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/unit/${NAME}")
    target_link_libraries(${NAME} PRIVATE ICEL FiNeR::FiNeR)
    target_include_directories(${NAME} PRIVATE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
    if(USE_TECIO)
        target_link_libraries(${NAME} PRIVATE -pthread)
    endif()
    if(USE_OPENMP)
        target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_Fortran)
    endif()
    if(USE_MPI)
        target_compile_definitions(${NAME} PRIVATE USE_MPI)
        target_link_libraries(${NAME} PRIVATE MPI::MPI_Fortran)
    endif()
    add_test(NAME ${NAME} COMMAND ${NAME} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/unit/${NAME}")
    set_tests_properties(${NAME} PROPERTIES LABELS "${ARGN}" TIMEOUT 120 ENVIRONMENT "${_ICE_TEST_ENV}")
endfunction()

# The property-table reader: header grammar, node and column checks, Tmin > 1
ice_add_unit_test(PropertiesUnit unit/test_properties.f90 "fast;unit;implementation")

# The solidification functions: the energy against IGLOO's hSolid, the recalescence, melting at T-melt, the
# nucleated fraction's events and majority threshold
ice_add_unit_test(SolidificationUnit unit/test_solidification.f90 "fast;unit;implementation")

# Rank count must not change the answer either. Needs an MPI build, so it is
# registered only when there is one.
if(USE_MPI)
    ice_add_fast_script(MPIEquivConnection fast/mpi-equiv "fast;implementation;needs-mpi")
    set_tests_properties(MPIEquivConnection PROPERTIES
        ENVIRONMENT "${_ICE_TEST_ENV};ICE_FAST_CASE=Doisneau/IG-split4")

    ice_add_fast_script(MPIEquivChimera fast/mpi-equiv "fast;implementation;chimera;needs-mpi")
    set_tests_properties(MPIEquivChimera PROPERTIES
        ENVIRONMENT "${_ICE_TEST_ENV};ICE_FAST_CASE=Doisneau/IG-chimera")

    # The same multigrid case again, this time also across ranks: each level carries
    # its own halo schedule, so a level switch has to move the exchange with it.
    ice_add_fast_script(MPIEquivMultigrid fast/multigrid "fast;implementation;multigrid;needs-mpi")
    set_tests_properties(MPIEquivMultigrid PROPERTIES
        ENVIRONMENT "${_ICE_TEST_ENV};ICE_MG_RANKS=2")

    # The two-block solidifying stream of SolidAdvection (leg O6) on two ranks, against the serial runs
    ice_add_verification(SolidAdvectionMPI O-solid-advection "fast;verification;phase-change;needs-mpi")
    set_tests_properties(SolidAdvectionMPI PROPERTIES
        ENVIRONMENT "${_ICE_TEST_ENV};ICE_SOLID_LEGS=O6;ICE_MPI_RANKS=2"
        RESOURCE_LOCK O-solid-advection-work)
endif()
