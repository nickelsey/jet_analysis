## build library for efficiency classes - initially, only
## run 4 as a stand in, until we get our embedding finalized
## for run 14

ADD_LIBRARY(Efficiency SHARED run4_eff.cc run14_eff.cc run7_eff.cc run7_eff_old.cc)
target_link_libraries(Efficiency ${ROOT_LIBRARIES})

#########################################
## install targets
INSTALL(TARGETS Efficiency DESTINATION lib/)
