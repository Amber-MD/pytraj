#ifndef INC_REFERENCEACTION_H
#define INC_REFERENCEACTION_H
#include "Trajin_Single.h"
#include "ReferenceFrame.h"
/// Class that can be used by Actions to hold reference structure/trajectory.
class ReferenceAction {
  public:
    ReferenceAction() {}
    /// Set up selected ref coordinates base on given frame and ref mask.
    inline void SetRefStructure(Frame const&, bool, bool);
    /// Initialize
    int InitRef(bool, bool, bool, bool, std::string const&, ReferenceFrame const&,
                Topology*, std::string const&, ArgList&, const char*);
    /// Setup
    int SetupRef(Topology const&, int, const char*);
    /// Perform necessary reference action based on mode
    inline void ActionRef(Frame const&, bool, bool);

    bool Previous()             const { return previous_;           } 
    const char* RefModeString() const { return modeString_.c_str(); }
    Frame const& RefFrame()     const { return refFrame_;           }
    Frame const& SelectedRef()  const { return selectedRef_;        }
    Vec3 const& RefTrans()      const { return refTrans_;           }
  private:
    enum RefModeType { UNKNOWN_REF=0, FIRST, REFFRAME, REFTRAJ };
    /// Set up ref mask for given topology. Allocate space for selected ref atoms.
    int SetRefMask(Topology const&, const char*);

    RefModeType refmode_;
    Frame refFrame_;         ///< Reference frame
    Frame selectedRef_;      ///< Atoms from reference frame selected by mask.
    AtomMask refMask_;       ///< Atoms to use from reference
    Vec3 refTrans_;          ///< If fitting, translation from origin to original ref center
    Trajin_Single refTraj_;  ///< Reference trajectory.
    bool previous_;          ///< True if current reference is previous frame (only RMSD now)
    std::string modeString_; ///< Information on current reference mode.
};
/// ----- INLINE FUNCTIONS -----------------------------------------------------
// ReferenceAction::SetRefStructure()
void ReferenceAction::SetRefStructure(Frame const& frameIn, bool fitIn, bool useMassIn)
{
  refFrame_ = frameIn;
  selectedRef_.SetCoordinates( refFrame_, refMask_ );
  if (fitIn)
    refTrans_ = selectedRef_.CenterOnOrigin( useMassIn );
}
// ReferenceAction::ActionRef()
void ReferenceAction::ActionRef(Frame const& frameIn, bool fitIn, bool useMassIn)
{
  if (refmode_ == FIRST) {
    SetRefStructure( frameIn, fitIn, useMassIn );
    refmode_ = REFFRAME;
  } else if (refmode_ == REFTRAJ) {
    refTraj_.GetNextFrame( refFrame_ );
    selectedRef_.SetCoordinates(refFrame_, refMask_);
    if (fitIn)
      refTrans_ = selectedRef_.CenterOnOrigin(useMassIn);
  }
}
#endif
