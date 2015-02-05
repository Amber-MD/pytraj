#ifndef INC_CPPTRAJSTATE_H
#define INC_CPPTRAJSTATE_H
#include "TrajinList.h"
#include "TrajoutList.h"
#include "TopologyList.h"
#include "DataSetList.h"
#include "DataFileList.h"
#include "ActionList.h"
#include "AnalysisList.h"
/// Hold all cpptraj state data
class CpptrajState {
  public:
    CpptrajState() : activeRef_(0), debug_(0), refidx_(0),
                     showProgress_(true), exitOnError_(true) {}
    // TODO: Change to &
    TopologyList* PFL()      { return &parmFileList_; }
    DataSetList* DSL()       { return &DSL_;          }
    DataFileList* DFL()      { return &DFL_;          }
    void SetNoExitOnError()  { exitOnError_ = false;  }
    void SetNoProgress()     { showProgress_ = false; }
    int Debug()        const { return debug_;         }
    bool ExitOnError() const { return exitOnError_;   }
    bool EmptyState()  const { return (actionList_.Empty() && 
                                       analysisList_.Empty() &&
                                       trajoutList_.Empty()); }
    void SetActionSilence(bool b)  { actionList_.SetSilent(b); }
    void SetActiveReference(DataSet_Coords_REF* rp) { activeRef_ = rp; }
    int AddTrajin( ArgList&, bool );
    int AddTrajin( std::string const& );
    int RunAnalyses();
    TrajinList const& InputTrajList() const { return trajinList_; }
    inline int AddTrajout( ArgList const& );
    inline int AddTrajout( std::string const& );
    int AddReference( std::string const&, ArgList const& );
    inline int AddReference( std::string const& );
    inline int AddAction( DispatchObject::DispatchAllocatorType, ArgList& );
    inline int AddAnalysis( DispatchObject::DispatchAllocatorType, ArgList& );
    static int WorldSize();
    static std::string PrintListKeys();
    int ListAll(ArgList&) const;
    int SetListDebug(ArgList&);
    int ClearList(ArgList&);
    int RemoveDataSet(ArgList&);
    int TrajLength( std::string const&, std::vector<std::string> const&);
    int Run();
    /// Write all DataFiles
    void MasterDataFileWrite();
  private:
    /// Types of lists
    enum ListType {
      L_ACTION = 0, L_TRAJIN, L_REF, L_TRAJOUT, L_PARM, L_ANALYSIS,
      L_DATAFILE, L_DATASET, N_LISTS
    };
    /// Hold list keyword.
    struct ListKeyType {
      ListType Type_;
      const char* Key_;
    };
    static ListKeyType ListKeys[];
    std::vector<bool> ListsFromArg(ArgList&, bool) const;
    /// \return active reference structure or empty frame if no reference.
    Frame ActiveReference() const;
    void ReferenceInfo() const;

    int RunNormal();
    int RunEnsemble();
    // -------------------------------------------
    /// List of parameter files 
    TopologyList parmFileList_;
     /// List of generated data sets
    DataSetList DSL_;
    /// List of datafiles that data sets will be written to
    DataFileList DFL_;
    /// List of input trajectory files
    TrajinList trajinList_;
    // -------------------------------------------
    /// List of actions to be performed each frame
    ActionList actionList_;
    /// List of output trajectory files 
    TrajoutList trajoutList_;
    // -------------------------------------------
    /// List of analyses to be performed on datasets
    AnalysisList analysisList_;
    
    /// Active reference structure for distance-based masks etc.
    DataSet_Coords_REF* activeRef_;
    /// State debug level
    int debug_;
    /// Internal reference index for numbering refs in DataSetList.
    int refidx_;
    /// Display Progress bar during run
    bool showProgress_;
    /// If true cpptraj will exit if errors are encountered instead of trying to continue
    bool exitOnError_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// CpptrajState::AddTrajout()
int CpptrajState::AddTrajout( ArgList const& argIn ) {
  return trajoutList_.AddTrajout( argIn, parmFileList_ );
}
int CpptrajState::AddTrajout( std::string const& fname ) {
  return AddTrajout( ArgList(fname) );
}
int CpptrajState::AddReference( std::string const& fname ) {
  return AddReference( fname, ArgList() );
}
// CpptrajState::AddAction()
int CpptrajState::AddAction( DispatchObject::DispatchAllocatorType Alloc, ArgList& argIn ) {
  argIn.MarkArg(0);
  return actionList_.AddAction( Alloc, argIn, &parmFileList_, &DSL_, &DFL_ );
}
// CpptrajState::AddAnalysis()
int CpptrajState::AddAnalysis( DispatchObject::DispatchAllocatorType Alloc, ArgList& argIn ) {
  argIn.MarkArg(0);
  return analysisList_.AddAnalysis( Alloc, argIn, &parmFileList_, &DSL_, &DFL_ );
}
#endif
