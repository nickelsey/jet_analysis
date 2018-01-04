## analysis scripts for the dijet imbalance measurement

## base dijet imbalance script
SET( BASIC_DIJET_IMBALANCE_BASE_SRC basic_dijet_imbalance.cc )
ADD_EXECUTABLE( basic_dijet_imbalance ${BASIC_DIJET_IMBALANCE_BASE_SRC} )
TARGET_LINK_LIBRARIES( basic_dijet_imbalance ${ROOT_LIBRARIES} ${STARPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} )
SET_TARGET_PROPERTIES( basic_dijet_imbalance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance/ )

## differential dijet imbalance analysis
## using the DijetWorker
SET( DIJET_IMBALANCE_BASE_SRC differential_aj.cc )
ADD_EXECUTABLE( dijet_imbalance ${DIJET_IMBALANCE_BASE_SRC} )
TARGET_LINK_LIBRARIES( dijet_imbalance ${ROOT_LIBRARIES} ${STARPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker Efficiency )
SET_TARGET_PROPERTIES( dijet_imbalance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance/ )


#########################################
## install targets
INSTALL( TARGETS basic_dijet_imbalance
         DESTINATION bin/dijet_imbalance/ )
INSTALL( TARGETS dijet_imbalance
         DESTINATION bin/dijet_imbalance/ )
