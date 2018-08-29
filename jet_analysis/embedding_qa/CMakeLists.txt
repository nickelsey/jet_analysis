## analysis scripts for data quality analysis

## generate efficiency curves from embedding QA
set(EFF_CURVE_SRCS make_eff_curves.cc)
add_executable(make_eff_curves ${EFF_CURVE_SRCS})
target_link_libraries(make_eff_curves ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency Centrality ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})
set_target_properties(make_eff_curves PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

#########################################
## install targets
install(TARGETS make_eff_curves
        DESTINATION bin/embed_qa)
