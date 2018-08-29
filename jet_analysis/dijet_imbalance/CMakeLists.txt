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
set(DIJET_IMBALANCE_BASE_SRC_AU differential_aj_auau.cc)
add_executable(dijet_imbalance_auau ${DIJET_IMBALANCE_BASE_SRC_AU})
target_link_libraries(dijet_imbalance_auau ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker Efficiency Centrality)
set_target_properties(dijet_imbalance_auau PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## differential dijet imbalance analysis for PP embedded in year 14
## using the DijetWorker
set(DIJET_IMBALANCE_BASE_SRC_PP differential_aj_pp.cc)
add_executable(dijet_imbalance_pp ${DIJET_IMBALANCE_BASE_SRC_PP})
target_link_libraries(dijet_imbalance_pp ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker Efficiency Centrality)
set_target_properties(dijet_imbalance_pp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## differential dijet imbalance analysis for PP embedded in year 7
## using the DijetWorker
set(DIJET_IMBALANCE_BASE_SRC_PP_Y7 differential_aj_pp_y7.cc)
add_executable(dijet_imbalance_pp_y7 ${DIJET_IMBALANCE_BASE_SRC_PP_Y7})
target_link_libraries(dijet_imbalance_pp_y7 ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker Efficiency Centrality)
set_target_properties(dijet_imbalance_pp_y7 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)


## given trees generated from the differential
## dijet imbalance measurement, produce Aj plots
## and various QA histograms
set(PRINT_AUAU_PP_RESULTS_SRC print_auau_pp_results.cc)
add_executable(print_auau_pp_results ${PRINT_AUAU_PP_RESULTS_SRC})
target_link_libraries(print_auau_pp_results ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(print_auau_pp_results PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## given trees generated from the differential
## dijet imbalance measurement, produce Aj plots
## and various QA histograms, including systematic errors
set(PRINT_AUAU_PP_SYS_SRC print_results_with_systematics.cc)
add_executable(print_results_with_systematics ${PRINT_AUAU_PP_SYS_SRC})
target_link_libraries(print_results_with_systematics ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(print_results_with_systematics PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## given trees generated from the differential
## dijet imbalance measurement, produce Aj plots
## and various QA histograms, including systematic errors
set(PRINT_AUAU_PP_SYS_QM_SRC print_results_with_systematics_qm.cc)
add_executable(print_results_with_systematics_qm ${PRINT_AUAU_PP_SYS_QM_SRC})
target_link_libraries(print_results_with_systematics_qm ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(print_results_with_systematics_qm PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

## test macro
set(COMPARE_WORKER_SRCS compare_worker_output.cc)
add_executable(compare_worker_output ${COMPARE_WORKER_SRCS})
target_link_libraries(compare_worker_output ${ROOT_LIBRARIES} ${Boost_LIBRARIES}  Efficiency Centrality ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})
set_target_properties(compare_worker_output PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/dijet_imbalance)

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
install(TARGETS dijet_imbalance_pp_y7
        DESTINATION bin/dijet_imbalance)
install(TARGETS print_auau_pp_results
        DESTINATION bin/dijet_imbalance)
install(TARGETS print_results_with_systematics
        DESTINATION bin/dijet_imbalance)
install(TARGETS print_results_with_systematics_qm
        DESTINATION bin/dijet_imbalance)
install(TARGETS compare_worker_output
        DESTINATION bin/dijet_imbalance)
