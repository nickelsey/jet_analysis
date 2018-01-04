## here we build a library out of the DijetMatrix class
## and its dependencies

ADD_LIBRARY( DijetWorker SHARED jet_def.cc match_def.hh dijet_definition.hh dijet_matrix.cc dijet_worker.cc )
TARGET_LINK_LIBRARIES( DijetWorker ${FASTJET_LIBRARIES} )

#########################################
## install targets
INSTALL( TARGETS DijetWorker DESTINATION lib/ )
