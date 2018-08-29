## here we build a library out of the DijetMatrix class
## and its dependencies

add_library(DijetWorker SHARED dijet_matrix/jet_def.cc dijet_matrix/match_def.cc dijet_matrix/dijet_definition.cc dijet_matrix/dijet_matrix.cc dijet_worker.cc off_axis_worker.cc)
target_link_libraries(DijetWorker ${FASTJET_LIBRARIES})

#########################################
## install targets
install(TARGETS DijetWorker DESTINATION lib)
