## analysis scripts for the dijet imbalance measurement

## base dijet imbalance script
set(BASIC_DIJET_IMBALANCE_BASE_SRC basic_dijet_imbalance.cc)
add_executable(basic_dijet_imbalance ${BASIC_DIJET_IMBALANCE_BASE_SRC})
target_link_libraries(basic_dijet_imbalance ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(basic_dijet_imbalance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## differential dijet imbalance analysis
## using the DijetWorker
set(DIJET_IMBALANCE_BASE_SRC differential_aj.cc)
add_executable(dijet_imbalance ${DIJET_IMBALANCE_BASE_SRC})
target_link_libraries(dijet_imbalance ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker Efficiency)
set_target_properties(dijet_imbalance PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## differential dijet imbalance analysis for AuAu
## using the DijetWorker
set(DIJET_IMBALANCE_BASE_SRC differential_aj_auau.cc)
add_executable(dijet_imbalance_auau ${DIJET_IMBALANCE_BASE_SRC})
target_link_libraries(dijet_imbalance_auau ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker Efficiency Centrality)
set_target_properties(dijet_imbalance_auau PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## differential dijet imbalance analysis for PP
## using the DijetWorker
set(DIJET_IMBALANCE_BASE_SRC differential_aj_pp.cc)
add_executable(dijet_imbalance_pp ${DIJET_IMBALANCE_BASE_SRC})
target_link_libraries(dijet_imbalance_pp ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker Efficiency Centrality)
set_target_properties(dijet_imbalance_pp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## given trees generated from the differential
## dijet imbalance measurement, produce Aj plots
## and various QA histograms
set(PRINT_AUAU_PP_RESULTS_SRC print_auau_pp_results.cc)
add_executable(print_auau_pp_results ${PRINT_AUAU_PP_RESULTS_SRC})
target_link_libraries(print_auau_pp_results ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(print_auau_pp_results PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

#########################################
## install targets
install(TARGETS basic_dijet_imbalance
         DESTINATION bin/dijet_imbalance)
install(TARGETS dijet_imbalance
         DESTINATION bin/dijet_imbalance)
install(TARGETS dijet_imbalance_auau
        DESTINATION bin/dijet_imbalance)
install(TARGETS dijet_imbalance_pp
        DESTINATION bin/dijet_imbalance)
install(TARGETS print_auau_pp_results
        DESTINATION bin/dijet_imbalance)
