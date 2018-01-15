## analysis scripts for various monte carlo routines

## testing AJ with different efficiencies
set(PYTHIA_AJ_EFFIC_SRC pythia_aj_efficiency.cc)
add_executable(pythia_aj_differential_efficiency ${PYTHIA_AJ_EFFIC_SRC})
target_link_libraries(pythia_aj_differential_efficiency ${ROOT_LIBRARIES} ${PYTHIA8_LIBRARIES} ${FASTJET_LIBRARIES} ${Boost_LIBRARIES} DijetWorker)
set_target_properties(pythia_aj_differential_efficiency PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/monte_carlo)

## print results from pythia_aj_differential_efficiency
set(PRINT_PYTHIA_AJ_EFFIC_SRC print_pythia_eff_results.cc)
add_executable(print_pythia_aj_efficiency_results ${PRINT_PYTHIA_AJ_EFFIC_SRC})
target_link_libraries(print_pythia_aj_efficiency_results ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(print_pythia_aj_efficiency_results PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/monte_carlo)

#########################################
## install targets
install(TARGETS pythia_aj_differential_efficiency
        DESTINATION bin/monte_carlo)
install(TARGETS print_pythia_aj_efficiency_results
        DESTINATION bin/monte_carlo)
