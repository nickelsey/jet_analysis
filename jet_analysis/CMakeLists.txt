## add the worker & utility functions to the include path

INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}/jet_analysis )

## add subdirectories to the project

ADD_SUBDIRECTORY( dijet_worker )
ADD_SUBDIRECTORY( dijet_imbalance )
ADD_SUBDIRECTORY( test )
ADD_SUBDIRECTORY( submit )
ADD_SUBDIRECTORY( efficiency )

