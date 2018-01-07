## here we build a library out of the DijetMatrix class
## and its dependencies

add_library(DijetWorker SHARED jet_def.cc match_def.hh dijet_definition.hh dijet_matrix.cc dijet_worker.cc)
target_link_libraries(DijetWorker ${FASTJET_LIBRARIES})

#########################################
## install targets
install(TARGETS DijetWorker DESTINATION lib)
