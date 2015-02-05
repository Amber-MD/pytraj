#ifndef INC_ANALYSIS_RMSAVGCORR_H
#define INC_ANALYSIS_RMSAVGCORR_H
#include "Analysis.h"
#include "DataSet_Coords.h"
// Class: Analysis_RmsAvgCorr
/// Calculate rmsd using running avg structures
class Analysis_RmsAvgCorr: public Analysis {
  public:
    Analysis_RmsAvgCorr();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_RmsAvgCorr(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    AtomMask tgtMask_;
    std::string separateName_;
    DataSet_Coords* coords_;
    DataSet* Ct_;
    DataSet* Csd_;
    Frame refFrame_;
    int maxwindow_;
    int lagOffset_;
    bool useMass_;
    bool useFirst_; ///< If true, use first running-avgd frame as reference.
};
#endif  
