#
# Distributed Linear Algebra Interface (DLAI)
#
# Copyright (c) 2018-2021, ETH Zurich
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause
#

# *** If (DLAI_INTERNAL_DOCUMENTATION) is ON !!! ***

# check if Doxygen is installed
find_package(Doxygen)

macro(transform_to_doxygen_list list)
  string(REPLACE ";" " " ${list} "${${list}}")
endmacro()

if (DOXYGEN_FOUND)
	message("Doxygen installed on system
	")
	
	set(DLAI_DOXY_PROJECT_NAME	"DLAI")
	
	# Doxygen source inputs
	set(DOXYGEN_INPUTS
		# Description
		${CMAKE_CURRENT_SOURCE_DIR}/Mainpage.dox
		# Installation documentation
		${CMAKE_CURRENT_SOURCE_DIR}/Start.dox
		# Examples & Tests
		${CMAKE_CURRENT_SOURCE_DIR}/Examples.dox
		# Inline documentation files
		${CMAKE_CURRENT_SOURCE_DIR}/../include/
		#${CMAKE_CURRENT_SOURCE_DIR}/../src/
		#${CMAKE_CURRENT_SOURCE_DIR}/../doc/
		# License text
		${CMAKE_CURRENT_SOURCE_DIR}/License.dox
	)
	transform_to_doxygen_list(DOXYGEN_INPUTS)

	set(DOXYGEN_EXAMPLE_PATH
		# doc path
		${CMAKE_CURRENT_SOURCE_DIR}/
		# test path
		${CMAKE_CURRENT_SOURCE_DIR}/../test/src/
	)
	transform_to_doxygen_list(DOXYGEN_EXAMPLE_PATH)
	
    # Configuration input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/../doc/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM
	)
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
