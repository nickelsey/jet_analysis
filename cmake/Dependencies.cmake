## list of project dependencies

########################################
## Set up ROOT environment - package is
## Built into ROOT6, I didn't check for ROOT5
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
message(STATUS "Found ROOT")

########################################
## Setup fastjet environment
find_package(FastJet REQUIRED)

#######################################
## Setup pythia environment
find_package(Pythia8 REQUIRED)

#######################################
## Setup boost environment
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS filesystem)
include_directories(${Boost_INCLUDE_DIRS})

#######################################
## Setup TStarJetPico environment
find_package(TStarJetPico REQUIRED)
