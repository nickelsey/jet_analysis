## build library for centrality classes - initially, only
## run 14

ADD_LIBRARY(Centrality SHARED centrality_run14.cc)

#########################################
## install targets
INSTALL(TARGETS Centrality DESTINATION lib/)

