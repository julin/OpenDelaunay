cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(ENABLE_OPENMP 1)
set(MODULE_NAME TestOpenDelaunay)

project (${MODULE_NAME})
 
set(DIR_PROJ_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(DIR_3P ${DIR_PROJ_ROOT}/3rd)
set(DIR_OUTPUT_ROOT ${DIR_PROJ_ROOT}/lib)
include(cmake/setting.cmake)

include(${DIR_PROJ_ROOT}/cmake/common.cmake)
add_subdirectory(3rd/googletest)

set(MODULE_SOURCES 
    src/delaunay3d.cxx
    src/sort_hilbert.cxx
    inc/delaunay3d.h
    test/test_config.h
)

include_directories(./)
include_directories(${DIR_PROJ_ROOT}/inc)
include_directories(${DIR_PROJ_ROOT}/3rd/googletest/googletest/include)

add_boost()
add_definitions(-D_ENABLE_UNITEST)

list(APPEND EXTERNAL_LIBRARIES gtest gtest_main)
add_executable(${MODULE_NAME} ${MODULE_SOURCES})

if(MSVC)
  set_target_properties(${MODULE_NAME} PROPERTIES COMPILE_FLAGS "/bigobj")
else()
  set_target_properties(${MODULE_NAME} PROPERTIES COMPILE_FLAGS "-std=c++11")
endif()

target_link_libraries(${MODULE_NAME} ${EXTERNAL_LIBRARIES} debug  "${EXTERNAL_LIBRARIES_DEBUG}" optimized  "${EXTERNAL_LIBRARIES_OPT}")


