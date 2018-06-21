// compare_auau.cc

#include "jet_analysis/util/common.hh"

DEFINE_string(test, "hello", "help message");

int main(int argc, char* argv[]) {
  
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  
  LOG(INFO) << "here?: " << FLAGS_test;
  
  
  return 0;
}
