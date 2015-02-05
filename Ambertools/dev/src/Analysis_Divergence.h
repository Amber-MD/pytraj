#ifndef INC_ANALYSIS_DIVERGENCE_H
#define INC_ANALYSIS_DIVERGENCE_H
#include "Analysis.h"
class Analysis_Divergence : public Analysis {
  public:
    Analysis_Divergence();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Divergence(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    /// Normalize sum over set to 1.0
    std::vector<double> NormalizeSet(DataSet const&, unsigned int) const;
    DataSet* ds1_;
    DataSet* ds2_;
};
#endif
