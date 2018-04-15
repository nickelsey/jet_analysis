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

#########################################
## install targets
install(TARGETS check_refmult comprehensive_data_quality
        DESTINATION bin/data_qa)
install(TARGETS print_data_quality generate_runid_list
        DESTINATION bin/data_qa)
install(TARGETS print_pp_compare generate_y14_efficiency
        DESTINATION bin/data_qa)
