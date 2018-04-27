# Usage: add_miniapp(miniapp_name [MPIRANKS <np>] [SRC <file.cpp>] [DEFS <defs>][LIBS <libs>])
# If EXE_NAME is not specified  ${miniapp_name}.cpp is used.
# If SRC is not specified  ${EXEC_NAME}.cpp is used.
# DEFS are the extra definitions needed to compile the test.
function(add_miniapp name)
  cmake_parse_arguments("arg" "" "EXE_NAME;SRC" "DEFS;LIBS" ${ARGN} )

  if(DEFINED arg_EXE_NAME)
    set(exe ${arg_EXE_NAME})
  else()
    set(exe ${name})
  endif()

  if(DEFINED arg_SRC)
    set(source ${arg_SRC})
  else()
    set(source ${exe}.cpp)
  endif()

  add_executable(${name} ${source})
  target_include_directories(${name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/test/include>)
  if(DEFINED arg_DEFS)
    target_compile_definitions(${name} PRIVATE ${arg_DEFS})
  endif()
  target_link_libraries(${name} gtest dla_interface ${arg_LIBS})
  set_target_properties(${name} PROPERTIES OUTPUT_NAME ${exe})

endfunction()
