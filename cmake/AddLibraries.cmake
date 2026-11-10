# ORION, as CMakeLists.txt required it with semver_dependency, and FiNeR (the submodule lib/FiNeR, with FiNeR's
# dependencies), each checked against the release ICE requires, so that a clone of another release stops the configure
# before anything is built. Included only when ICE is the top-level project: a parent project gives ICE its libraries.

include(${CMAKE_CURRENT_LIST_DIR}/SubmodulePins.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/FiNeRDependencies.cmake)

# ORION: a release selected by SemVerDeps satisfies the requirement; an external clone (local) only warns there,
# so ICE checks its version, the newest vX.Y.Z tag of its HEAD as SemVerDeps found it
get_property(_orion_source GLOBAL PROPERTY SEMVER_ORION_SOURCE)
if (_orion_source STREQUAL "local")
  get_property(_orion_version GLOBAL PROPERTY SEMVER_ORION_VERSION)
  set(_ok FALSE)
  if (_orion_version)
    _semver_satisfies("${_orion_version}" _ok ${ICE_ORION_VERSION})
  else()
    set(_orion_version "of no release (no vX.Y.Z tag in its history)")
  endif()
  if (NOT _ok)
    _semver_range("${ICE_ORION_VERSION}" _min _max)
    message(${ICE_PIN_CHECK} "ORION at ${FETCHCONTENT_SOURCE_DIR_ORION} is ${_orion_version}, but ICE requires "
                    "ORION ${ICE_ORION_VERSION} (${_min} up to ${_max} excluded). Give a clone of such a release, or "
                    "build the release ICE selects (without FETCHCONTENT_SOURCE_DIR_ORION), or configure with "
                    "-DICE_PIN_CHECK=WARNING to build this ORION on purpose.")
  endif()
endif()

# FiNeR: its dependencies before its CMakeLists adds them. Those of the submodule are fetched; an external
# clone is built as it is, its dependencies only checked. Its release is in its VERSION file (v2.3.0): its
# project() still declares 2.0.3
FetchContent_Declare(FiNeR SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/FiNeR")
if (FETCHCONTENT_SOURCE_DIR_FINER)
  set(ICE_FINER_SOURCE_DIR "${FETCHCONTENT_SOURCE_DIR_FINER}")
  ensure_finer_dependencies("${ICE_FINER_SOURCE_DIR}" FALSE)
else()
  set(ICE_FINER_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/FiNeR")
  ensure_submodule_initialized("${ICE_FINER_SOURCE_DIR}")
  ice_check_submodule_pin("${ICE_FINER_SOURCE_DIR}")
  ensure_finer_dependencies("${ICE_FINER_SOURCE_DIR}" TRUE)
endif()
set(_finer_version "")
if (EXISTS "${ICE_FINER_SOURCE_DIR}/VERSION")
  file(STRINGS "${ICE_FINER_SOURCE_DIR}/VERSION" _finer_version LIMIT_COUNT 1)
  string(REGEX REPLACE "^[ \t]*v" "" _finer_version "${_finer_version}")
endif()
string(REGEX REPLACE "^v" "" _finer_pinned "${ICE_FINER_TAG}")
string(CONCAT _hint "Give FETCHCONTENT_SOURCE_DIR_FINER a clone of FiNeR ${ICE_FINER_TAG}, or build the submodule "
                    "lib/FiNeR (without it); if lib/FiNeR was bumped on purpose, set ICE_FINER_TAG in "
                    "CMakeLists.txt to its release")
ice_check_version(FiNeR "${ICE_FINER_SOURCE_DIR}" "${_finer_version}" "${_finer_pinned}" "${_hint}")

FetchContent_MakeAvailable(ORION FiNeR)
