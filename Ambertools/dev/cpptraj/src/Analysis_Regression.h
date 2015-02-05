#ifndef INC_ANALYSIS_REGRESSION_H
#define INC_ANALYSIS_REGRESSION_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_Regression : public Analysis {
  public:
    Analysis_Regression();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Regression(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    Array1D output_dsets_;
    CpptrajFile outfile_;
};
#endif
