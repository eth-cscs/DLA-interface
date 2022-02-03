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

message("<FindDLAF.cmake>")

set(DLAI_WITH_DLAF_INTERNAL OFF)
	unset(DLAF_FOUND CACHE)

if (DLAF_ROOT)
	set(DLAF_DIR PATH ${DLAF_ROOT} "DLAF_DIR directory")
	message(STATUS "DLA-Future path: ${DLAF_DIR}")
else()
	message(STATUS "To enable DLAF must set DLAF_ROOT variable!")
endif()
  
if (DLAF_DIR)
	message(STATUS "DLAF_DIR = ${DLAF_DIT}")	
endif()

message("</FindDLAF.cmake>")