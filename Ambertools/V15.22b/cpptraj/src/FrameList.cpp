// FrameList
#include "FrameList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
FrameList::FrameList() : activeRefNum_(0), debug_(0) {}

// DESTRUCTOR
FrameList::~FrameList() {
  Clear();
}

/** Clear the FrameList. */
void FrameList::Clear() {
  // NOTE: Because ref frames are currently passed around by pointers they
  //       are deleted here instead of in ReferenceFrame; this avoids double
  //       frees when e.g. a ReferenceFrame is destroyed in an Action.
  // TODO: ReferenceFrame should be passed by const reference. 
  for (std::vector<ReferenceFrame>::iterator ref = frames_.begin();
                                             ref != frames_.end(); ++ref)
    ref->ClearRef();
  frames_.clear();
  activeRefNum_ = 0;
}

void FrameList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0)
    mprintf("ReferenceList debug level set to %i\n", debug_);
}

// -----------------------------------------------------------------------------
// FrameList::ActiveReference()
/** Return the address of the frame pointed to by activeRefNum_.
  */
Frame FrameList::ActiveReference() const {
  if (frames_.empty()) return Frame();
  return frames_[activeRefNum_].Coord();
}

// FrameList::SetActiveRef()
/** Set the given frame list number as the active reference.
  */
int FrameList::SetActiveRef(int numIn) {
  if (numIn < 0 || numIn >= (int)frames_.size()) {
    mprintf("Warning: FrameList::SetActiveRef: Ref # %i out of bounds.\n",numIn);
    return 1;
  }
  activeRefNum_ = numIn;
  return 0;
}
// -----------------------------------------------------------------------------
/** Add Frame from the given trajectory file to the FrameList. Store trajectory 
  * name and frame number that this frame came from in frameNames and frameNums 
  * respectively. Store the associated parm in FrameParm. 
  */
int FrameList::AddRefFrame(ArgList& argIn, TopologyList const& topListIn) {
  ReferenceFrame refFrame;

  // Get filename, associated parmtop, and mask.
  std::string fname = argIn.GetStringNext();
  Topology* parmIn = topListIn.GetParm( argIn );
  std::string mask  = argIn.GetMaskNext();
  // 'average' keyword is deprecated
  if ( argIn.hasKey("average") ) {
    mprinterr("Error: 'average' for reference is deprecated. Please use\n"
              "Error:   the 'average' action to create averaged coordinates.\n");
    return 1;
  }
  // Load reference frame.
  if (refFrame.LoadRef( fname, argIn, parmIn, mask, debug_ ))
    return 1;
  // Check and warn if filename/reftag currently in use
  // TODO: Check for base filename?
  for (std::vector<ReferenceFrame>::const_iterator rf = frames_.begin();
                                                   rf != frames_.end(); ++rf)
  {
    if ( rf->FrameName().Full() == refFrame.FrameName().Full() )
      mprintf("Warning: Reference with name '%s' already exists!\n",
              refFrame.FrameName().full());
    if ( !refFrame.Tag().empty() && rf->Tag() == refFrame.Tag() )
      mprintf("Error: Reference with tag [%s] already exists!\n",
              refFrame.Tag().c_str());
  }
  // Add reference frame to list.
  frames_.push_back( refFrame );

  return 0;
}

const char* FrameList::RefArgs = "reference | ref <name> | refindex <#>";

// FrameList::GetFrameFromArgs()
/** \return ReferenceFrame based on args in argIn.
  * The keywords in order of precedence are:
  *   - 'ref <name>'  : Get reference frame by name/tag.
  *   - 'reference'   : First reference frame in list.
  *   - 'refindex <#>': Reference frame at position.
  */
ReferenceFrame FrameList::GetFrameFromArgs(ArgList& argIn) const {
  int refindex;
  // By name/tag
  std::string refname = argIn.GetStringKey("ref");
  if (!refname.empty()) {
    ReferenceFrame rf = GetFrameByName( refname );
    if (rf.empty()) {
      mprinterr("Error: Could not get reference with name %s\n", refname.c_str());
      return ReferenceFrame(-1);
    }
    return rf; 
  }
  // First defined reference
  if (argIn.hasKey("reference")) {
    if (frames_.empty()) {
      mprinterr("Error: No reference frames defined.\n");
      return ReferenceFrame(-1);
    }
    return frames_[0];
  }
  // By index
  refindex = argIn.getKeyInt("refindex", -1);
  if (refindex != -1) {
    if (refindex < 0 || refindex >= (int)frames_.size()) {
      mprinterr("Error: reference index %i is out of bounds.\n", refindex);
      return ReferenceFrame(-1);
    }
    return frames_[refindex];
  }
  // No frame desired, return empty.
  return ReferenceFrame();
}

// FrameList::GetFrameByName()
ReferenceFrame FrameList::GetFrameByName(std::string const& refName) const {
  for (std::vector<ReferenceFrame>::const_iterator rf = frames_.begin();
                                                   rf != frames_.end(); ++rf)
  { // OK if name matches tag, full path, or base file name.
    if ( rf->Tag()==refName || rf->FrameName().MatchFullOrBase( refName ) )
      return *rf;
  }
  return ReferenceFrame();
}

// FrameList::List()
/** Print a list of trajectory names that frames have been taken from.
  */
void FrameList::List() const {
  if (!frames_.empty()) {
    mprintf("\nREFERENCE FRAMES (%zu total):\n", frames_.size());
    for (std::vector<ReferenceFrame>::const_iterator rf = frames_.begin();
                                                     rf != frames_.end(); ++rf)
    {
      mprintf("    %u:", rf - frames_.begin());
      rf->RefInfo();
    }
    mprintf("\tActive reference frame for masks is %i\n",activeRefNum_);
  }
}
