## analysis scripts for data quality analysis

## checks refmult distribution against
## hand calculated refmult
set(CHECK_REFMULT_SRCS check_refmult.cc)
add_executable(check_refmult ${CHECK_REFMULT_SRCS})
target_link_libraries(check_refmult ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(check_refmult PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## full data quality for event level, tracks level
## and calorimeter level QA
set(COMPREHENSIVE_DQ_SRCS comprehensive_data_quality.cc)
add_executable(comprehensive_data_quality ${COMPREHENSIVE_DQ_SRCS})
target_link_libraries(comprehensive_data_quality ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(comprehensive_data_quality PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

#########################################
## install targets
install(TARGETS check_refmult comprehensive_data_quality
        DESTINATION bin/data_qa)
