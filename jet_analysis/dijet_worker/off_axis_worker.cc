#include "jet_analysis/dijet_worker/off_axis_worker.hh"

OffAxisWorker::OffAxisWorker() {

}

OffAxisWorker::~OffAxisWorker() {

}

std::unordered_map<std::string, std::unique_ptr<OffAxisOutput>>& OffAxisWorker::Run(DijetWorker& worker) {
    
    // clear our output
    off_axis_result.clear();
    
    // get the input container
    auto& dijets = worker.Dijets();
    
    for (auto& dijet : dijets) {
        
    }
  
  return off_axis_result;
}
