# Sets DLA-Future variables.
# - Must set DLAF_ROOT variable
#
# DLAF_DIR contains the directory in which DLA-Future has Config file.
# DLAI_WITH_DLAF is set to ON if DLAF is available, OFF otherwise.

function(dla_find_dlaf)
	set(DLAI_WITH_DLAF_INTERNAL OFF)
	unset(DLAF_FOUND CACHE)

	if (DLAF_ROOT)
		setoption(DLAF_DIR PATH ${DLAF_ROOT}/lib64/cmake "DLAF_DIR directory")
	else()
		message(STATUS "To enable DLAF must set DLAFROOT variable!")
	endif()
  
	if (DLAF_DIR)
		find_package(DLAF)
		
		if(DLAF_FOUND)
			link_libraries(DLAF::DLAF)
			set(DLAI_WITH_DLAF_INTERNAL ON)
			set(DLAI_WITH_DLAF ${DLAI_WITH_DLAF_INTERNAL} CACHE BOOL "DLAF is available (autogenerated)" FORCE) 
			message(STATUS "Found DLA-Future: ${DLAFROOT}")
		else()
			message(STATUS "DLA-Future: NOT FOUND!")
		endif()
		
		if (DLAI_WITH_DLAF)
			add_compile_definitions(DLAI_WITH_DLAF)
			add_compile_definitions(DLAF_FUNCTION_NAME=__PRETTY_FUNCTION__)
			message(STATUS "DLA-Future: Found DLAF and set DLAI_WITH_DLAF!")
		endif()
	endif()
  
endfunction()
