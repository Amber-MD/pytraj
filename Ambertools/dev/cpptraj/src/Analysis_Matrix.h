#ifndef INC_ANALYSIS_MATRIX_H
#define INC_ANALYSIS_MATRIX_H
#include "Analysis.h"
#include "DataSet_2D.h"
#include "DataSet_Modes.h"
class Analysis_Matrix : public Analysis {
  public:
    Analysis_Matrix();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Matrix(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    int NMWizOutput() const;

    DataSet_2D* matrix_;
    DataSet_Modes* modes_;
    std::string outthermo_;
    double thermo_temp_;
    int nevec_;
    bool thermopt_;
    bool reduce_;
    bool eigenvaluesOnly_;
    bool nmwizopt_;
    int nmwizvecs_;
    std::string nmwizfile_;
    Topology nmwizParm_;
};
#endif
