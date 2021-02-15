if (DEFINED HWLOC_ROOT OR DEFINED ENV{HWLOC_ROOT})
  find_path(HWLOC_INCLUDE_DIRS hwloc.h REQUIRED)
  find_library(HWLOC_LIBRARIES hwloc REQUIRED)
else()
  find_package(PkgConfig REQUIRED QUIET)
  pkg_search_module(HWLOC hwloc IMPORTED_TARGET GLOBAL)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HWLOC DEFAULT_MSG
  HWLOC_INCLUDE_DIRS
  HWLOC_LIBRARIES
)

if (HWLOC_FOUND)
  if (NOT TARGET HWLOC::HWLOC)
    add_library(HWLOC::HWLOC INTERFACE IMPORTED GLOBAL)
  endif()

  target_include_directories(HWLOC::HWLOC INTERFACE ${HWLOC_INCLUDE_DIRS})
  target_link_libraries(HWLOC::HWLOC INTERFACE ${HWLOC_LIBRARIES})
endif()
