#ifndef OFF_AXIS_WORKER_HH
#define OFF_AXIS_WORKER_HH

#include "jet_analysis/dijet_worker/dijet_worker.hh"

struct OffAxisOutput {
    
};

class OffAxisWorker {
    public:
        OffAxisWorker();
        ~OffAxisWorker();

        std::unordered_map<std::string, std::unique_ptr<OffAxisOutput>>& Run(DijetWorker& worker);

        std::unordered_map<std::string, std::unique_ptr<OffAxisOutput>>& OffAxisResult() {return off_axis_result;}

    private:
        std::unordered_map<std::string, std::unique_ptr<OffAxisOutput>> off_axis_result;

};

#endif // OFF_AXIS_WORKER_HH
