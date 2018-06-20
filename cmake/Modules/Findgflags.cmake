# - Try to find GFLAGS
#
# The following variables are optionally searched for defaults
#  GFLAGS_ROOT_DIR:            Base directory where all GFLAGS components are found
#
# The following are set after configuration is done:
#  GFLAGS_FOUND
#  GFLAGS_INCLUDE_DIRS
#  GFLAGS_LIBRARIES
#  GFLAGS_LIBRARY_DIRS

include(FindPackageHandleStandardArgs)

set(GFLAGS_ROOT_DIR "" CACHE PATH "Folder contains Gflags")

# We are testing only a couple of files in the include directories
if(NOT WIN32)
    find_path(GFLAGS_INCLUDE_DIR gflags/gflags.h
        PATHS ${GFLAGS_ROOT_DIR})
endif()

if(MSVC)
    find_package(gflags NO_MODULE HINTS ${GFLAGS_ROOT_DIR})
    set(GFLAGS_LIBRARY ${gflags_LIBRARIES})
else()
    message(STATUS "GOT HERE: " ${GFLAGS_ROOT_DIR})
    find_library(GFLAGS_LIBRARY gflags
    PATHS ${GFLAGS_ROOT_DIR}
    PATH_SUFFIXES lib lib64)
    message(STATUS "LIBRARY: " ${GFLAGS_LIBRARY})
endif()

find_package_handle_standard_args(gflags DEFAULT_MSG GFLAGS_INCLUDE_DIR GFLAGS_LIBRARY)


if(GFLAGS_FOUND)
    set(GFLAGS_INCLUDE_DIRS ${GFLAGS_INCLUDE_DIR})
    set(GFLAGS_LIBRARIES ${GFLAGS_LIBRARY})
    message(STATUS "Found gflags  (include: ${GFLAGS_INCLUDE_DIR}, library: ${GFLAGS_LIBRARY})")
    mark_as_advanced(GFLAGS_LIBRARY_DEBUG GFLAGS_LIBRARY_RELEASE
                     GFLAGS_LIBRARY GFLAGS_INCLUDE_DIR GFLAGS_ROOT_DIR)
else(GFLAGS_FOUND)
    unset(GFLAGS_INCLUDE_DIRS CACHE)
    unset(GFLAGS_INCLUDE_DIR CACHE)
    unset(GFLAGS_LIBRARY CACHE)
    unset(GFLAGS_LIBRARIES CACHE)
    unset(GFLAGS_ROOT_DIR CACHE)
    unset(GFLAGS_LIBRARY_DIR CACHE)
endif()
