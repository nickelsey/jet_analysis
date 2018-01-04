###########################################
## root CMake file for jet analysis code ##
##                  WSU                  ##
###########################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.6)

PROJECT( jet_analysis )

## Set CMake search paths for user defined modules
## These are in the project's source dir
SET( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )

## don't let cmake install to /usr/local since we don't have root
## on the grid
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "default install path" FORCE )
endif()

## the grid has an older version of GCC - we have to explicitly
## turn on support for c++11
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11" )

########################################
## Set up ROOT environment - package is
## Built into ROOT6, I didn't check for ROOT5
LIST( APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} )
FIND_PACKAGE( ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net )
INCLUDE( ${ROOT_USE_FILE} )
MESSAGE( STATUS "Found ROOT" )

########################################
## Setup fastjet environment
FIND_PACKAGE( FastJet REQUIRED )

#######################################
## Setup pythia environment
FIND_PACKAGE( Pythia8 REQUIRED )

#######################################
## Setup boost environment
SET( Boost_USE_STATIC_LIBS OFF )
SET( Boost_USE_MULTITHREADED ON )
SET( Boost_USE_STATIC_RUNTIME OFF )
FIND_PACKAGE( Boost REQUIRED COMPONENTS filesystem )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )

#######################################
## Setup TStarJetPico environment
FIND_PACKAGE( StarPico REQUIRED )

#######################################
## enable testing
ENABLE_TESTING()

#######################################
## Including subdirectories
ADD_SUBDIRECTORY( jet_analysis/ )

#######################################
## install source files in install
## directory
INSTALL( DIRECTORY jet_analysis DESTINATION . 
         FILES_MATCHING PATTERN "CMakeLists.txt" EXCLUDE
         PATTERN "._*" EXCLUDE
         PATTERN "submit" EXCLUDE
         PATTERN "*")
