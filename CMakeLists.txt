###########################################
## root CMake file for jet analysis code ##
##                  WSU                  ##
###########################################

cmake_minimum_required(VERSION 3.6)

project(jet_analysis)

## don't let cmake install to /usr/local since we don't have root
## on the grid
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX
       "${CMAKE_BINARY_DIR}/install"
       CACHE PATH "default install path" FORCE)
endif()

## the grid has an older version of GCC - we have to explicitly
## turn on support for c++11
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")

#######################################
## lookup project dependencies
## add module path to cmake
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

#######################################
## dependencies
include(cmake/Dependencies.cmake)

#######################################
## system checks
include(cmake/SysCheck.cmake)

#######################################
## utility functions
include(cmake/Utils.cmake)

#######################################
## forcing proper include lookup paths - must take
## configured files in build dir before empty files
## in the source dir
include_directories(BEFORE ${PROJECT_SOURCE_DIR})
include_directories(BEFORE ${PROJECT_BINARY_DIR})

#######################################
## enable testing
enable_testing()

#######################################
## Including subdirectories
add_subdirectory(jet_analysis/)

#######################################
## install source files in install
## directory
install(DIRECTORY jet_analysis DESTINATION .
        FILES_MATCHING PATTERN "CMakeLists.txt" EXCLUDE
        PATTERN "._*" EXCLUDE
        PATTERN "*.in.*" EXCLUDE
        PATTERN "submit" EXCLUDE
        PATTERN "*")

## add an uninstall option, since cmake does not provide
## one by default
if(NOT TARGET uninstall)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
