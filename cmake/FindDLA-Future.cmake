#
# CMake recipes
#
# Copyright (c) 2018-2021, ETH Zurich
# BSD 3-Clause License. All rights reserved.
#
# author: Aleksander GRM (aleksander.grm@fpp.uni-lj.si)
#

# Find DLA-Future library
#
# ScaLAPACK depends on MPI and LAPACK, so it depends on other modules for their respective
# targets LAPACK::LAPACK and MPI::MPI_<LANG>. In particular, for the latter one, this module checks
# which language is enabled in the project and it adds all needed dependencies.
#
# Users can manually specify next variables (even by setting them empty to force use of
# the compiler implicit linking) to control which implementation they want to use:
#   SCALAPACK_LIBRARY
#       ;-list of {lib name, lib filepath, -Llibrary_folder}
#
# This module sets the following variables:
#   SCALAPACK_FOUND - set to true if a library implementing the SCALAPACK interface is found
#
# If ScaLAPACK symbols got found, it creates target SCALAPACK::SCALAPACK

cmake_minimum_required(VERSION 3.12)


if (DLAF_ROOT)
#	message(STATUS "DLA-Future path: ${DLAF_ROOT}")
	find_package(DLAF)
else()
	message(STATUS "Found DLA-Future: FALSE")
	message(STATUS "To enable DLAF must set DLAF_ROOT variable to point to DLA-Future install path! (spack location -i dla-future)")
endif()

if (DLAF_FOUND)
	message(STATUS "Found DLA-Future: TRUE (found version hpx)")
	link_libraries(DLAF::DLAF)
	set(DLAI_WITH_DLAF ON)
	add_compile_definitions(DLAI_WITH_DLAF)
	add_compile_definitions(DLAF_FUNCTION_NAME=__PRETTY_FUNCTION__)
else()
	message(STATUS "DLA-Future NOT FOUND!")
endif()