# google test framework
#set(gtest_inc ${PROJECT_SOURCE_DIR}/test/gtest)
add_library(gtest gtest/src/gtest-all.cpp gtest/src/gtest_main.cpp)
target_include_directories(gtest PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gtest>)

# Usage: add_test(src_name_without_extension [libs])
function(add_unit_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} gtest ${argv})
  add_test(NAME ${name} COMMAND "$<TARGET_FILE:${name}>")
endfunction()

add_subdirectory(unit)
