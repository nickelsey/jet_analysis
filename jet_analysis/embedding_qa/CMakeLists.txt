## analysis scripts for data quality analysis

## generate efficiency curves from embedding QA
set(EFF_CURVE_SRCS make_eff_curves.cc)
add_executable(make_eff_curves ${EFF_CURVE_SRCS})
target_link_libraries(make_eff_curves ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency Centrality ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})
set_target_properties(make_eff_curves PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## generate efficiency curves from embedding QA scaled by DCA differences
set(EFF_CURVE_SCALED_SRCS make_eff_curves_scaled.cc)
add_executable(make_eff_curves_scaled ${EFF_CURVE_SCALED_SRCS})
target_link_libraries(make_eff_curves_scaled ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency Centrality ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})
set_target_properties(make_eff_curves_scaled PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)

## generate efficiency curves from embedding QA
set(Y7_PT_SRCS make_y7_pt_spectra.cc)
add_executable(make_y7_pt_spectra ${Y7_PT_SRCS})
target_link_libraries(make_y7_pt_spectra ${ROOT_LIBRARIES} ${TSTARJETPICO_LIBRARIES} ${Boost_LIBRARIES} Efficiency Centrality ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES})
set_target_properties(make_y7_pt_spectra PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/data_qa)



#########################################
## install targets
install(TARGETS make_eff_curves make_eff_curves_scaled
        DESTINATION bin/embed_qa)
install(TARGETS make_y7_pt_spectra
        DESTINATION bin/embed_qa)
