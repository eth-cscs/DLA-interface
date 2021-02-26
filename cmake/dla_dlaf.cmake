# Sets DLA-Future variables.
# - Must set DLAFROOT variable
#
# DLAF_DIR contains the directory in which DLA-Future has Config file.
# DLA_HAVE_DLAF is set to ON if DLAF is available, OFF otherwise.

function(dla_find_dlaf)
	set(DLA_HAVE_DLAF_INTERNAL OFF)
	unset(DLAF_FOUND CACHE)

	if (DLAFROOT)
		setoption(DLAF_DIR PATH ${DLAFROOT}/lib/cmake "DLAF_DIR directory")
	else()
		message(STATUS "To enable DLAF must set DLAFROOT variable!")
	endif()
  
	if (DLAF_DIR)
		find_package(DLAF)
		
		if(DLAF_FOUND)
			link_libraries(DLAF::DLAF)
			set(DLA_HAVE_DLAF_INTERNAL ON)
			set(DLA_HAVE_DLAF ${DLA_HAVE_DLAF_INTERNAL} CACHE BOOL "DLAF is available (autogenerated)" FORCE) 
			message(STATUS "Found DLA-Future: ${DLAFROOT}")
		else()
			message(STATUS "DLA-Future: NOT FOUND!")
		endif()
		
		if (DLA_HAVE_DLAF)
			add_compile_definitions(DLA_HAVE_DLAF)
			add_compile_definitions(DLAF_FUNCTION_NAME=__PRETTY_FUNCTION__)
			message(STATUS "DLA-Future: Found DLAF and set DLA_HAVE_DLAF!")
		endif()
	endif()
  
endfunction()
