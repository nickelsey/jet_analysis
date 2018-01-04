########################################
## Set up TStarJetPico library environment
## STARPICODIR or STARPICOPATH have to be set as an environment variable
#  STARPICO_FOUND - System has TStarJetPico
#  STARPICO_INCLUDE_DIRS - The TStarJetPico include directories
#  STARPICO_LIBRARIES - The libraries of the TStarJetPico framework


# - Try to find TStarJetPico
# Once done this will define

find_path( STARPICO_INCLUDE_DIR TStarJetPicoReader.h
           HINTS $ENV{STARPICOPATH} $ENV{STARPICODIR} )

find_library( STARPICO_LIBRARY NAMES TStarJetPico
              HINTS $ENV{STARPICODIR} $ENV{STARPICOPATH} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set STARPICO_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args( StarPico
                                   STARPICO_INCLUDE_DIR STARPICO_LIBRARY )

mark_as_advanced( STARPICO_LIBRARY STARPICO_INCLUDE_DIR )

set(STARPICO_LIBRARIES ${STARPICO_LIBRARY} )
set(STARPICO_INCLUDE_DIRS ${STARPICO_INCLUDE_DIR} )

## Add TStarJetPico include path to cmake
INCLUDE_DIRECTORIES ( ${STARPICO_INCLUDE_DIRS} )
LINK_DIRECTORIES ( ${STARPICO_INCLUDE_DIRS} )
