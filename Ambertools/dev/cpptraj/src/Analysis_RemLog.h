#ifndef INC_ANALYSIS_REMLOG_H
#define INC_ANALYSIS_REMLOG_H
#include "Analysis.h"
#include "DataSet_RemLog.h"
class Analysis_RemLog : public Analysis {
  public:
    Analysis_RemLog();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_RemLog(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    enum ModeType { NONE = 0, CRDIDX, REPIDX };
    bool calculateStats_;
    bool calculateLifetimes_;
    bool printIndividualTrips_;
    DataSet_RemLog* remlog_;
    ModeType mode_;
    std::vector<DataSet*> outputDsets_;
    CpptrajFile statsout_;
    CpptrajFile reptime_;
    std::string lifetimesName_;
    std::string acceptout_;
    int calcRepFracSlope_;
    CpptrajFile repFracSlope_;
};
#endif
