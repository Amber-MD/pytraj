#ifndef INC_ANALYSIS_CROSSCORR_H
#define INC_ANALYSIS_CROSSCORR_H
#include "Analysis.h"
class Analysis_CrossCorr : public Analysis {
  public:
    Analysis_CrossCorr();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_CrossCorr(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();

  private:
    DataFile* outfile_;
    DataSetList dsets_;
    DataSet* matrix_;
};
#endif
