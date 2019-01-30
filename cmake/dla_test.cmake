# Usage: add_test(test_name [MPIRANKS <np>] [SRC <file.cpp>] [DEFS <defs>][LIBS <libs>])
# If SRC is not specified  ${test_name}.cpp is used.
# DEFS are the extra definitions needed to compile the test.
# If MPIRANKS is defined the test is linked with the MPI library and runned with
# ${TEST_RUNNER} ${RUNNER_NP_OPT} ${np} <exe>
function(add_unit_test name)
  cmake_parse_arguments("arg" "" "MPIRANKS" "SRC;DEFS;LIBS" ${ARGN} )

  if(DEFINED arg_SRC)
    set(source ${arg_SRC})
  else()
    set(source ${name}.cpp)
  endif()

  add_executable(${name} ${source})
  target_include_directories(${name} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/test/include>)
  if(DEFINED arg_DEFS)
    target_compile_definitions(${name} PRIVATE ${arg_DEFS})
  endif()
  target_link_libraries(${name} gtest dla_interface ${arg_LIBS})

  if(DEFINED arg_MPIRANKS)
    if(NOT arg_MPIRANKS GREATER 0)
      message(FATAL_ERROR "Wrong MPIRANKS number ${arg_MPIRANKS}")
    endif()
    target_link_libraries(${name} ${MPI_CXX_LIBRARIES})
    add_test(NAME ${name} COMMAND ${TEST_RUNNER} ${TEST_RUNNER_NP_OPT} ${arg_MPIRANKS} "$<TARGET_FILE:${name}>")
  else()
    if(DLA_ALL_TESTS_USE_RUNNER)
      add_test(NAME ${name} COMMAND ${TEST_RUNNER} ${TEST_RUNNER_NP_OPT} 1 "$<TARGET_FILE:${name}>")
    else()
      add_test(NAME ${name} COMMAND "$<TARGET_FILE:${name}>")
    endif()
  endif()
endfunction()
