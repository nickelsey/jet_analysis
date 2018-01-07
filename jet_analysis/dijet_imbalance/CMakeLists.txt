## analysis scripts for the dijet imbalance measurement

## base dijet imbalance script
set(BASIC_DIJET_IMBALANCE_BASE_SRC basic_dijet_imbalance.cc)
add_executable(basic_dijet_imbalance ${BASIC_DIJET_IMBALANCE_BASE_SRC})
target_link_libraries(basic_dijet_imbalance ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(basic_dijet_imbalance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance/)

## differential dijet imbalance analysis
## using the DijetWorker
set(DIJET_IMBALANCE_BASE_SRC differential_aj.cc)
add_executable(dijet_imbalance ${DIJET_IMBALANCE_BASE_SRC})
target_link_libraries(dijet_imbalance ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker Efficiency)
set_target_properties(dijet_imbalance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance/)


#########################################
## install targets
install(TARGETS basic_dijet_imbalance
         DESTINATION bin/dijet_imbalance)
install(TARGETS dijet_imbalance
         DESTINATION bin/dijet_imbalance)
