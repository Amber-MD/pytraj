#ifndef INC_ANALYSISLIST_H
#define INC_ANALYSISLIST_H
#include "Analysis.h"
/// Hold all analyses to be performed.
class AnalysisList {
  public:
    AnalysisList();
    ~AnalysisList();
    void Clear(); 
    void SetDebug(int);
    int Debug() const { return debug_; }
    int AddAnalysis(DispatchObject::DispatchAllocatorType, ArgList&,
                    TopologyList*, DataSetList*, DataFileList*);
    int DoAnalyses();
    void List() const;
    bool Empty() const { return analysisList_.empty(); }
  private:
    /// Analysis setup status
    enum AnalysisStatusType { NO_SETUP = 0, SETUP, INACTIVE };
    typedef std::vector<Analysis*> aListType;
    /// List of analyses
    aListType analysisList_;
    /// List of analysis commands
    std::vector<std::string> analysisCmd_;
    /// List of analysis statuses
    std::vector<AnalysisStatusType> analysisStatus_;
    /// Default debug level for analyses
    int debug_;
};
#endif
