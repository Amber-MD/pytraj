#ifndef INC_ANALYSIS_FFT_H
#define INC_ANALYSIS_FFT_H
#include "Analysis.h"
#include "Array1D.h"
/// Calculate FFT of dataset(s)
class Analysis_FFT : public Analysis {
  public:
    Analysis_FFT();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_FFT(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    Array1D input_dsets_;
    Array1D output_dsets_;
    double dt_; ///< Sampling interval (timestep)
};
#endif
