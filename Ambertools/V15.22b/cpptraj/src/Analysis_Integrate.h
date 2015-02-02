#ifndef INC_ANALYSIS_INTEGRATE_H
#define INC_ANALYSIS_INTEGRATE_H
#include "Analysis.h"
#include "Array1D.h"
#include "DataSet_Mesh.h"
class Analysis_Integrate : public Analysis {
  public:
    Analysis_Integrate();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Integrate(); }
    static void Help();
  
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataFile* outfile_; // FIXME: May not need to be class var
    Array1D input_dsets_;
    std::vector<DataSet_Mesh*> output_dsets_;
};
#endif
