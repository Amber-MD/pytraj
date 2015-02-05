#ifndef INC_ANALYSIS_AMDBIAS_H
#define INC_ANALYSIS_AMDBIAS_H
#include "Analysis.h"
class Analysis_AmdBias : public Analysis {
  public:
    Analysis_AmdBias();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_AmdBias(); }
    static void Help();
  
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSet* ds1_;
    double Ethresh_;
    double alpha_;
    DataSet* bias_;
};
#endif
