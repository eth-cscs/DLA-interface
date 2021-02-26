
include(dla_utils)
include(CheckFunctionExists)

function(dla_find_hwloc)
  unset(HWLOC_FOUND CACHE)
  setoption(HWLOC_ROOT PATH "${HWLOCROOT}" "HWLOC root directory")
  
  if (HWLOC_ROOT)
    set(CMAKE_PREFIX_PATH ${HWLOC_ROOT})
    find_package(PkgConfig REQUIRED)
    pkg_search_module(HWLOC REQUIRED hwloc)
    unset(CMAKE_PREFIX_PATH)
    set(DLA_HWLOC_LIBS "-L${HWLOC_LIBDIR};${HWLOC_LIBRARIES}" CACHE STRING "Libraries for HWLOC" FORCE)
    add_definitions(-DDLA_HAVE_HWLOC)
  else()
    message(WARNING "Configurations without HWLOC may have bad performance due to bad thread binding. Please refers to the DLA interface documentation.")
  endif()

endfunction()
