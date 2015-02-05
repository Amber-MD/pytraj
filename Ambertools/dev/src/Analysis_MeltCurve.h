#ifndef INC_ANALYSIS_MELTCURVE_H
#define INC_ANALYSIS_MELTCURVE_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_MeltCurve : public Analysis {
  public:
    Analysis_MeltCurve() : mcurve_(0), cut_(0.0) {}
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_MeltCurve(); }
    static void Help();
  
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    DataSet* mcurve_;
    double cut_;
};
#endif
