#ifndef INC_ANALYSIS_SPLINE_H
#define INC_ANALYSIS_SPLINE_H
#include "Analysis.h"
#include "Array1D.h"
#include "DataSet_Mesh.h"
class Analysis_Spline : public Analysis {
  public:
    Analysis_Spline();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Spline(); }
    static void Help();
  
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    DataFile* outfile_; // FIXME: May not need to be class var
    Array1D input_dsets_;
    std::vector<DataSet_Mesh*> output_dsets_;
    int meshsize_;       ///< Mesh size will be DataSet.Size * meshfactor
    double meshmin_;     ///< Min value of resulting mesh
    double meshmax_;     ///< Max value of resulting mesh
    double meshfactor_;  ///< Min value of input data sets
    bool useDefaultMin_; ///< If true use meshmin_ for all meshes.
    bool useDefaultMax_; ///< If true use meshmax_ for all meshes. 
};
#endif
