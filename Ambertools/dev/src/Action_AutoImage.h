#ifndef INC_ACTION_AUTOIMAGE_H
#define INC_ACTION_AUTOIMAGE_H
#include "Action.h"
class Action_AutoImage : public Action {
  public:
    Action_AutoImage();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_AutoImage(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    AtomMask anchorMask_; ///< Used to center anchor region.
    std::string anchor_;  ///< Mask expression for anchor region.
    std::string fixed_;   ///< Mask expression for fixed region.
    std::string mobile_;  ///< Mask expression for mobile region.
    bool origin_;         ///< If true imaging occurs w.r.t. coordinate origin.
    bool ortho_;          ///< If true imaging is orthogonal.
    bool usecom_;         ///< If true imaging of mobile region uses molecule center.
    bool truncoct_;
    bool useMass_;
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic_; ///< Determine whether triclinic code should be used.

    typedef std::vector<int> pairList;
    pairList anchorList_;
    pairList fixedList_;
    pairList mobileList_;

    pairList SetupAtomRanges(Topology*,std::string const&);
};
#endif
