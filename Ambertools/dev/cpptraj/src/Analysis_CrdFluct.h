#ifndef INC_ANALYSIS_CRDFLUCT_H
#define INC_ANALYSIS_CRDFLUCT_H
#include "Analysis.h"
#include "DataSet_Coords.h"
class Analysis_CrdFluct : public Analysis {
  public:
    Analysis_CrdFluct();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_CrdFluct(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataSet_Coords* coords_;
    AtomMask mask_;
    typedef std::vector<DataSet*> SetList;
    SetList outSets_;
    bool bfactor_;
    int windowSize_;

    void CalcBfactors( Frame, Frame, double, DataSet& );
};
#endif
