#ifndef INC_ANALYSIS_H
#define INC_ANALYSIS_H
#include "DispatchObject.h"
#include "ArgList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "TopologyList.h"
// Class: Analysis
/// The abstract base class that all other analyses inherit.
/** Analysis occurs after trajectories are read and data sets populated.
  * Analysis operates on those data sets.
  */
class Analysis : public DispatchObject {
  public:
    /// Enumerate potential return stats from Setup and Analyze.
    enum RetType { OK = 0, ERR };
    /// Destructor - virtual since this class is inherited
    virtual ~Analysis() {}
    /// Set up action
    virtual RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int) = 0;
    virtual RetType Analyze() = 0;
};
#endif
