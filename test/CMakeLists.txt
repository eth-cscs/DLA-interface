# google test framework
#set(gtest_inc ${PROJECT_SOURCE_DIR}/test/gtest)
add_library(gtest gtest/src/gtest-all.cpp)
target_include_directories(gtest PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gtest>)

add_library(gtest_main gtest/src/gtest_main.cpp)
target_include_directories(gtest_main PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gtest>)

include(dla_test)

add_subdirectory(unit)
