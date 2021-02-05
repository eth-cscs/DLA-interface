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
		find_package(DLAF CONFIG)
		
		if(DLAF_FOUND)
			set(DLAF_LINK_DIRECTORIES ${DLAFROOT}/lib CACHE STRING "Link paths for DLAF")
			set(DLA_DLAF_INCLUDE_DIRS ${DLAFROOT}/include CACHE STRING "Include paths for DLAF")
			set(DLA_DLAF_INCLUDE_DIRS ${DLAFROOT}/include CACHE STRING "Include paths for DLAF")
			set(DLA_DLAF_LIBRARIES DLAF::DLAF CACHE STRING "Libraries for DLAF")
			
			set(CMAKE_REQUIRED_INCLUDES ${DLA_DLAF_INCLUDE_DIRS})
			set(CMAKE_REQUIRED_LIBRARIES ${DLA_DLAF_LIBRARIES})
			set(DLA_HAVE_DLAF_INTERNAL ON)
			set(DLA_HAVE_DLAF ${DLA_HAVE_DLAF_INTERNAL} CACHE BOOL "DLAF is available (autogenerated)" FORCE) 
			message(STATUS "Found DLA-Future: ${DLAFROOT}")
			
			link_directories(${DLAF_LINK_DIRECTORIES})
			include_directories(${DLA_DLAF_INCLUDE_DIRS};${DLAFROOT}/../blaspp/include;${DLAFROOT}/../lapackpp/include;${DLAFROOT}/../hpx/include)
			
			message(DEBUG "Set LINK_DIRECTORIES: ${DLAF_LINK_DIRECTORIES}")
			message(DEBUG "Set INCLUDE_DIR: ${DLA_DLAF_INCLUDE_DIRS}")
			message(DEBUG "Set LIBRARIES: ${DLA_DLAF_LIBRARIES}")
		else()
			message(STATUS "DLA-Future: NOT FOUND!")
		endif()
		
		if (DLA_HAVE_DLAF)
			add_compile_definitions(DLA_HAVE_DLAF)
			add_compile_definitions(DLAF_FUNCTION_NAME=__PRETTY_FUNCTION__)
		endif()
	endif()
  
endfunction()
