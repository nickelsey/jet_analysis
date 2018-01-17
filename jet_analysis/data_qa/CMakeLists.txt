## analysis scripts for data quality analysis

## base dijet imbalance script
set(CHECK_REFMULT_SRCS check_refmult.cc)
add_executable(check_refmult ${CHECK_REFMULT_SRCS})
target_link_libraries(check_refmult ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(check_refmult PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

#########################################
## install targets
install(TARGETS check_refmult
        DESTINATION bin/data_qa)
