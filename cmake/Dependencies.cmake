## list of project dependencies

#######################################
## Set up ROOT environment - package is
## Built into ROOT6, I didn't check for ROOT5
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
message(STATUS "Found ROOT")

#######################################
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

#######################################
## setup glog and gflags
find_package(gflags REQUIRED)
find_package(glog REQUIRED)

#######################################
## testing resources
## googletest
set(TEMP_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF)
set(BUILD_GTEST ON)
set(INSTALL_GTEST OFF)
## gmock currently not used
set(BUILD_GMOCK OFF)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/googletest)
include_directories(${PROJECT_SOURCE_DIR}/third_party/googletest/googletest/include)

## benchmark
# We will not need to test benchmark lib itself.
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark testing.")
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "Disable benchmark install to avoid overwriting.")
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/benchmark)
include_directories(${PROJECT_SOURCE_DIR}/third_party/benchmark/include)

# restore the build shared libs option.
set(BUILD_SHARED_LIBS ${TEMP_BUILD_SHARED_LIBS})

