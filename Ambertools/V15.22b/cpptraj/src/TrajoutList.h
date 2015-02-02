#ifndef INC_TRAJOUTLIST_H
#define INC_TRAJOUTLIST_H
#include "Trajout.h"
#include "TopologyList.h"
// Class: TrajoutList
/// Hold trajectories for output
class TrajoutList {
  public:
    TrajoutList();
    ~TrajoutList();
    void Clear();
    void SetDebug(int);
    int AddEnsembleTrajout(ArgList const&, TopologyList const&, int);
    /// Add a traj file to the list with given access and associate with a parm
    int AddTrajout(ArgList const&, TopologyList const&);
    /// Call write for all trajectories
    int WriteTrajout(int, Topology*, Frame*);
    /// Call end for all trajectories
    void CloseTrajout();
    void List() const;
    bool Empty()     const { return trajout_.empty();     }
    // The definitions below are for ensemble processing.
    typedef std::vector<ArgList> ArgsArray;
    typedef std::vector<ArgList>::const_iterator ArgIt;
    ArgIt argbegin() const { return trajoutArgs_.begin(); }
    ArgIt argend()   const { return trajoutArgs_.end();   }
  private:
    int debug_;
    typedef std::vector<Trajout*> ListType;
    ListType trajout_;
    /// Array of trajout args for setting up ensemble trajout.
    ArgsArray trajoutArgs_;

    int AddTrajout(std::string const&, ArgList&, TopologyList const&,
                   TrajectoryFile::TrajFormatType);
};
#endif
