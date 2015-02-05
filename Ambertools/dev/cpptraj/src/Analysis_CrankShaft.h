#ifndef INC_ANALYSIS_CRANKSHAFT_H
#define INC_ANALYSIS_CRANKSHAFT_H
#include "Analysis.h"
#include "DataSet_1D.h"
class Analysis_CrankShaft : public Analysis {
  public: 
    Analysis_CrankShaft();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_CrankShaft(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    enum CStype { ANGLE=0, DISTANCE };
    static const char* CSstring[];
    static const char* distance_ss_2D[][6];
    static const char* torsion_ss_2D[][6];

    std::string filename_;    ///< Filename for frame vs bin number
    std::string resultsname_; ///< Filename for crankshaft results
    int debug_;
    int start_;
    int stop_;
    int offset_;
    CStype type_;
    DataSet_1D* scalar1_;
    DataSet_1D* scalar2_;
    std::string info_;
};
#endif
