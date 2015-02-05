#ifndef INC_TRAJINLIST_H
#define INC_TRAJINLIST_H
#include "Trajin.h"
#include "TopologyList.h"
// Class: TrajinList
/// Hold input trajectories
class TrajinList {
  public:
    enum TrajModeType { UNDEFINED, NORMAL, ENSEMBLE };
    TrajinList();
    ~TrajinList();
    void Clear();
    void SetDebug(int dIn) { debug_ = dIn; }
    /// Add a traj file to the list based on input from arg list
    int AddTrajin(std::string const&, ArgList&, TopologyList const&);
    /// Add trajectory ensemble
    int AddEnsemble(std::string const&, ArgList&, TopologyList const&);

    typedef std::vector<Trajin*> ListType;
    typedef ListType::const_iterator const_iterator;
    const_iterator begin() const { return trajin_.begin(); }
    const_iterator end()   const { return trajin_.end();   }
    bool empty()           const { return trajin_.empty(); }
    TrajModeType Mode()    const { return mode_;           }
    const Trajin* front()  const { return trajin_.front(); }
    int MaxFrames()        const { return maxframes_;      }
    void List() const;
  private:
    int AddInputTraj(std::string const&, Trajin*, ArgList, TopologyList const&);

    ListType trajin_;
    int debug_;
    int maxframes_;
    TrajModeType mode_; ///< Trajectory processing mode
    /// CRDIDXARG: Used when processing ensemble and sorting by CRDIDX
    std::string finalCrdIndicesArg_;
};
#endif

