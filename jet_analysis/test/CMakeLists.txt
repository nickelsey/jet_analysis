#############################################
## tests, including simple tests to make sure
## everything has linked properly

## tests if ROOT include directory and libraries are set properly
SET( ROOT_BASIC_TESTING_SRCS basic_root_test.cc )
ADD_EXECUTABLE( root_test ${ROOT_BASIC_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( root_test ${ROOT_LIBRARIES} )
SET_TARGET_PROPERTIES( root_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

## tests if FastJet include directory and libraries are set properly
SET( FASTJET_TESTING_SRCS basic_fastjet_test.cc )
ADD_EXECUTABLE( fastjet_test ${FASTJET_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( fastjet_test ${FASTJET_LIBRARIES} ${ROOT_LIBRARIES} )
SET_TARGET_PROPERTIES( fastjet_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

## tests if TStarJetPico include directory and libraries are set properly
SET(TS_TESTING_SRCS basic_ts_test.cc )
ADD_EXECUTABLE( ts_test ${TS_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( ts_test ${ROOT_LIBRARIES} ${STARPICO_LIBRARIES} )
SET_TARGET_PROPERTIES( ts_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( JET_DEF_TESTING_SRCS jet_def_test.cc )
ADD_EXECUTABLE( jet_def_test ${JET_DEF_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( jet_def_test DijetWorker)
SET_TARGET_PROPERTIES( jet_def_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( MATCH_DEF_TESTING_SRCS match_def_test.cc )
ADD_EXECUTABLE( match_def_test ${MATCH_DEF_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( match_def_test DijetWorker )
SET_TARGET_PROPERTIES( match_def_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( DIJET_DEFINITION_TESTING_SRCS dijet_definition_test.cc )
ADD_EXECUTABLE( dijet_def_test ${DIJET_DEFINITION_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( dijet_def_test DijetWorker )
SET_TARGET_PROPERTIES( dijet_def_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( DIJET_MATRIX_TESTING_SRCS dijet_matrix_test.cc )
ADD_EXECUTABLE( dijet_matrix_test ${DIJET_MATRIX_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( dijet_matrix_test DijetWorker )
SET_TARGET_PROPERTIES( dijet_matrix_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( DIJET_WORKER_TESTING_SRCS dijet_worker_test.cc )
ADD_EXECUTABLE( dijet_worker_test ${DIJET_WORKER_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( dijet_worker_test DijetWorker )
SET_TARGET_PROPERTIES( dijet_worker_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( ARG_HELPER_TESTING_SRCS arg_helper_test.cc )
ADD_EXECUTABLE( arg_helper_test ${ARG_HELPER_TESTING_SRCS} )
SET_TARGET_PROPERTIES( arg_helper_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( TRIGGER_LOOKUP_TESTING_SRCS trigger_lookup_test.cc)
ADD_EXECUTABLE( trig_lookup_test ${TRIGGER_LOOKUP_TESTING_SRCS} )
SET_TARGET_PROPERTIES( trig_lookup_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( SELECTOR_COMPARE_TESTING_SRCS selector_compare_test.cc )
ADD_EXECUTABLE( selector_compare_test ${SELECTOR_COMPARE_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( selector_compare_test ${FASTJET_LIBRARIES} )
SET_TARGET_PROPERTIES( selector_compare_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( DIJET_KEY_TESTING_SRCS dijet_key_test.cc )
ADD_EXECUTABLE( dijet_key_test ${DIJET_KEY_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( dijet_key_test DijetWorker )
SET_TARGET_PROPERTIES( dijet_key_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

## reader util testing has a few text files it depends on,
## which are copied into the binary directory
CONFIGURE_FILE( csv_test.csv ${CMAKE_BINARY_DIR}/bin/test/csv_test.csv )
CONFIGURE_FILE( reader_util_test.in.txt ${CMAKE_SOURCE_DIR}/jet_analysis/test/reader_util_test.txt )
CONFIGURE_FILE( reader_util_test.in.txt ${CMAKE_BINARY_DIR}/bin/test/reader_util_test.txt )
CONFIGURE_FILE( example_bad_run_list.txt ${CMAKE_BINARY_DIR}/bin/test/example_bad_run_list.txt )
CONFIGURE_FILE( example_bad_tower_list.txt ${CMAKE_BINARY_DIR}/bin/test/example_bad_tower_list.txt )

## and also must be configured itself
CONFIGURE_FILE( reader_util_test.in.cc ${CMAKE_SOURCE_DIR}/jet_analysis/test/reader_util_test.cc )

SET( READER_UTIL_TESTING_SRCS reader_util_test.cc )
ADD_EXECUTABLE( reader_util_test ${READER_UTIL_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( reader_util_test ${ROOT_LIBRARIES} ${STARPICO_LIBRARIES} ${Boost_LIBRARIES} )
SET_TARGET_PROPERTIES( reader_util_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )

SET( RUN4_EFF_TESTING_SRCS run4_eff_test.cc )
ADD_EXECUTABLE( run4_eff_test ${RUN4_EFF_TESTING_SRCS} )
TARGET_LINK_LIBRARIES( run4_eff_test ${ROOT_LIBRARIES} Efficiency )
SET_TARGET_PROPERTIES( run4_eff_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/test/ )


#########################################
## Add tests
ADD_TEST( NAME star_basic COMMAND ts_test  )
ADD_TEST( NAME fastjet_basic COMMAND fastjet_test )
ADD_TEST( NAME root_basic COMMAND root_test )
ADD_TEST( NAME arg_helper COMMAND arg_helper_test )
ADD_TEST( NAME trig_lookup COMMAND trig_lookup_test )
ADD_TEST( NAME jet_def COMMAND jet_def_test )
ADD_TEST( NAME match_def COMMAND match_def_test )
ADD_TEST( NAME dijet_definition COMMAND dijet_def_test )
ADD_TEST( NAME selector_compare COMMAND selector_compare_test )
ADD_TEST( NAME dijet_key COMMAND dijet_key_test )
ADD_TEST( NAME dijet_matrix COMMAND dijet_matrix_test )
ADD_TEST( NAME dijet_worker COMMAND dijet_worker_test )
ADD_TEST( NAME reader_util COMMAND reader_util_test )
ADD_TEST( NAME run4_eff COMMAND run4_eff_test )

#########################################
## install targets
INSTALL( TARGETS ts_test fastjet_test root_test
         DESTINATION bin/test/ )
INSTALL( TARGETS jet_def_test arg_helper_test
         DESTINATION bin/test/ )
INSTALL( TARGETS match_def_test dijet_def_test
         DESTINATION bin/test/ )
INSTALL( TARGETS selector_compare_test dijet_key_test
         DESTINATION bin/test/ )
INSTALL( TARGETS dijet_matrix_test dijet_worker_test
         DESTINATION bin/test/ )
INSTALL( TARGETS reader_util_test run4_eff_test
         DESTINATION bin/test/ )
