#include "Action_Unwrap.h"
#include "CpptrajStdio.h"
#include "ImageRoutines.h"

// CONSTRUCTOR
Action_Unwrap::Action_Unwrap() :
  imageMode_(Image::BYATOM),
  RefParm_(0),
  orthogonal_(false),
  center_(false)
{ }

void Action_Unwrap::Help() {
  mprintf("\t[center] [{bymol | byres | byatom}]\n"
          "\t[ %s ] [<mask>]\n", FrameList::RefArgs);
  mprintf("  Reverse of 'image'; unwrap coordinates in <mask> according\n"
          "  to a reference structure.\n");
}

// Action_Unwrap::Init()
Action::RetType Action_Unwrap::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get Keywords
  center_ = actionArgs.hasKey("center");
  if (actionArgs.hasKey("bymol"))
    imageMode_ = Image::BYMOL;
  else if (actionArgs.hasKey("byres"))
    imageMode_ = Image::BYRES;
  else if (actionArgs.hasKey("byatom")) {
    imageMode_ = Image::BYATOM;
    // Unwrapping to center by atom makes no sense
    if (center_) center_ = false;
  } else
    imageMode_ = Image::BYATOM;
  // Get reference
  ReferenceFrame REF = FL->GetFrameFromArgs( actionArgs );
  if (REF.error()) return Action::ERR;
  if (!REF.empty()) {
    RefFrame_ = REF.Coord();
    // Get reference parm for frame
    RefParm_ = (Topology*)(&REF.Parm()); // FIXME: 
  }

  // Get mask string
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    UNWRAP: By %s using mask '%s'", 
          Image::ModeString(imageMode_), mask_.MaskString());
  if (imageMode_ != Image::BYATOM) {
    if (center_)
      mprintf(" based on center of mass.");
    else
      mprintf(" based on first atom position.");
  }
  mprintf("\n");
  if ( !REF.empty())
    mprintf("\tReference is %s", REF.FrameName().base());
  else
    mprintf("\tReference is first frame.");
  mprintf("\n");

  return Action::OK;
}

// Action_Unwrap::Setup()
Action::RetType Action_Unwrap::Setup(Topology* currentParm, Topology** parmAddress) {
  // Ensure same number of atoms in current parm and ref parm
  if ( RefParm_!=0 ) {
    if ( currentParm->Natom() != RefParm_->Natom() ) {
      mprinterr("Error: unwrap: # atoms in reference parm %s is not\n", RefParm_->c_str());
      mprinterr("Error:         equal to # atoms in parm %s\n", currentParm->c_str());
      return Action::ERR;
    }
  }
  // Check box type
  if (currentParm->BoxType()==Box::NOBOX) {
    mprintf("Error: unwrap: Parm %s does not contain box information.\n",
            currentParm->c_str());
    return Action::ERR;
  }
  orthogonal_ = false;
  if (currentParm->BoxType()==Box::ORTHO)
    orthogonal_ = true;

  // Setup atom pairs to be unwrapped.
  imageList_ = Image::CreatePairList(*currentParm, imageMode_, mask_);
  if (imageList_.empty()) {
    mprintf("Warning: Mask '%s' selects no atoms for topology '%s'.\n",
            mask_.MaskString(), currentParm->c_str());
    return Action::ERR;
  }
  mprintf("\tNumber of %ss to be unwrapped is %zu based on mask '%s'\n",
          Image::ModeString(imageMode_), imageList_.size()/2, mask_.MaskString());

  // Use current parm as reference if not already set
  if (RefParm_ == 0)
    RefParm_ = currentParm;
  return Action::OK;
}

// Action_Unwrap::DoAction()
Action::RetType Action_Unwrap::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 ucell, recip;
  // Set reference structure if not already set
  if (RefFrame_.empty()) {
    RefFrame_ = *currentFrame;
    return Action::OK;
  }
 
  if (orthogonal_)
    Image::UnwrapOrtho( *currentFrame, RefFrame_, imageList_, center_, true );
  else {
    currentFrame->BoxCrd().ToRecip( ucell, recip );
    Image::UnwrapNonortho( *currentFrame, RefFrame_, imageList_, ucell, recip, center_, true );
  }

  return Action::OK;
} 
