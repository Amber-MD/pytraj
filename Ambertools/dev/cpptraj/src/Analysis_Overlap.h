#ifndef INC_ANALYSIS_OVERLAP_H
#define INC_ANALYSIS_OVERLAP_H
#include "Analysis.h"
class Analysis_Overlap : public Analysis {
  public:
    Analysis_Overlap();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Overlap(); }
    static void Help();
  
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSet* ds1_;
    DataSet* ds2_;
    bool useDeviation_;
};
#endif
