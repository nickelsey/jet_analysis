## analysis scripts for data quality analysis

## checks refmult distribution against
## hand calculated refmult
set(CHECK_REFMULT_SRCS check_refmult.cc)
add_executable(check_refmult ${CHECK_REFMULT_SRCS})
target_link_libraries(check_refmult ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(check_refmult PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## creates a list of all runIDs present in a data file
set(GENERATE_RUNID_LIST_SRCS generate_runid_list)
add_executable(generate_runid_list ${GENERATE_RUNID_LIST_SRCS})
target_link_libraries(generate_runid_list ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(generate_runid_list PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)


## full data quality for event level, tracks level
## and calorimeter level QA
set(COMPREHENSIVE_DQ_SRCS comprehensive_data_quality.cc)
add_executable(comprehensive_data_quality ${COMPREHENSIVE_DQ_SRCS})
target_link_libraries(comprehensive_data_quality ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency)
set_target_properties(comprehensive_data_quality PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## print routines to make comparisons between luminosities
## and run periods
set(DQ_PRINT_SRCS print_data_quality.cc)
add_executable(print_data_quality ${DQ_PRINT_SRCS})
target_link_libraries(print_data_quality ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(print_data_quality PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## print routines to make comparisons between pp data sets
set(DQ_PP_PRINT_SRCS print_compare_pp.cc)
add_executable(print_pp_compare ${DQ_PP_PRINT_SRCS})
target_link_libraries(print_pp_compare ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(print_pp_compare PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## generate efficiency curves
set(GEN_EFF_SRCS generate_y14_efficiency.cc)
add_executable(generate_y14_efficiency ${GEN_EFF_SRCS})
target_link_libraries(generate_y14_efficiency ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(generate_y14_efficiency PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## generate pt spectra and QA for different data sets
set(GEN_PT_SPEC_SRCS generate_pt_spectra.cc)
add_executable(generate_pt_spectra ${GEN_PT_SPEC_SRCS})
target_link_libraries(generate_pt_spectra ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency Centrality)
set_target_properties(generate_pt_spectra PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## compare pt spectra for y14 & y7
set(COMPARE_PT_SRCS compare_y14_y7_pt.cc)
add_executable(compare_y14_y7_pt ${COMPARE_PT_SRCS})
target_link_libraries(compare_y14_y7_pt ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency Centrality)
set_target_properties(compare_y14_y7_pt PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## compare different AuAu datasets
set(COMPARE_AUAU_SRCS compare_auau.cc)
add_executable(compare_auau ${COMPARE_AUAU_SRCS})
target_link_libraries(compare_auau ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} Efficiency Centrality)
set_target_properties(compare_auau PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## compare different AuAu datasets - print
set(COMPARE_AUAU_PRINT_SRCS compare_auau_print.cc)
add_executable(compare_auau_print ${COMPARE_AUAU_PRINT_SRCS})
target_link_libraries(compare_auau_print ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} Efficiency Centrality)
set_target_properties(compare_auau_print PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## create nprimary-nglobal cut for a given dataset
set(GLOBAL_PRIMARY_SRCS global_primary.cc)
add_executable(global_primary ${GLOBAL_PRIMARY_SRCS})
target_link_libraries(global_primary ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} Efficiency Centrality)
set_target_properties(global_primary PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## compare pt spectra for y14 & y7
set(NGLOBAL_PT_SRCS nglobal_cut_pt_spectra.cc)
add_executable(nglobal_cut_pt_spectra ${NGLOBAL_PT_SRCS})
target_link_libraries(nglobal_cut_pt_spectra ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency Centrality)
set_target_properties(nglobal_cut_pt_spectra PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## compare P16id and P17id production
set(P17_P16_SRCS p16id_p17id_compare.cc)
add_executable(p16id_p17id_compare ${P17_P16_SRCS})
target_link_libraries(p16id_p17id_compare ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency Centrality ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})
set_target_properties(p16id_p17id_compare PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)


#########################################
## install targets
install(TARGETS check_refmult comprehensive_data_quality
        DESTINATION bin/data_qa)
install(TARGETS print_data_quality generate_runid_list
        DESTINATION bin/data_qa)
install(TARGETS print_pp_compare generate_y14_efficiency
        DESTINATION bin/data_qa)
install(TARGETS generate_pt_spectra compare_y14_y7_pt
        DESTINATION bin/data_qa)
install(TARGETS compare_auau compare_auau_print
        DESTINATION bin/data_qa)
install(TARGETS global_primary nglobal_cut_pt_spectra
        DESTINATION bin/data_qa)
install(TARGETS p16id_p17id_compare
        DESTINATION bin/data_qa)
