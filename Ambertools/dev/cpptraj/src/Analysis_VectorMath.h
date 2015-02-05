#ifndef INC_ANALYSIS_VECTORMATH_H
#define INC_ANALYSIS_VECTORMATH_H
#include "Analysis.h"
#include "DataSet_Vector.h"
class Analysis_VectorMath : public Analysis {
  public:
    Analysis_VectorMath();
    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_VectorMath(); }
    static void Help();
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    enum ModeType {  DOTPRODUCT = 0, DOTANGLE, CROSSPRODUCT };
    static const char* ModeString[];

    int DotProduct();
    int CrossProduct();

    ModeType mode_;
    DataSet_Vector* vinfo1_;
    DataSet_Vector* vinfo2_;
    DataSet* DataOut_;       ///< Output data set
    bool norm_;
};
#endif
