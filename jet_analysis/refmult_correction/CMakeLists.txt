## analyses related to creating refmult corrections
## for STAR data

## generate corrections for refmult distribution
## as a function of luminosity & vz
set(LUMI_VZ_CORRECTION_SRCS lumi_vz_correction.cc)
add_executable(lumi_vz_correction ${LUMI_VZ_CORRECTION_SRCS})
target_link_libraries(lumi_vz_correction ${ROOT_LIBRARIES} ${Boost_LIBRARIES})
set_target_properties(lumi_vz_correction PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/refmult_correction)

#########################################
## install targets
install(TARGETS lumi_vz_correction
        DESTINATION bin/refmult_correction)
