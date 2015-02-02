// Action_Center 
#include "Action_Center.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Center::Action_Center() :
  centerMode_(Frame::BOXCTR),
  useMass_(false)
{ } 

void Action_Center::Help() {
  mprintf("\t<mask> [origin] [mass]\n"
          "\t [ %s [<refmask>]]\n  Center coordinates in <mask>.\n", FrameList::RefArgs);
}

// Action_Center::Init()
Action::RetType Action_Center::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  if (actionArgs.hasKey("origin"))
    centerMode_ = Frame::ORIGIN;
  else
    centerMode_ = Frame::BOXCTR;
  useMass_ = actionArgs.hasKey("mass");
  ReferenceFrame refFrm = FL->GetFrameFromArgs( actionArgs );
  if (refFrm.error()) return Action::ERR;

  // Get Masks
  Mask_.SetMaskString( actionArgs.GetMaskNext() );
  // Get reference mask if reference specified.
  AtomMask refMask;
  if (!refFrm.empty()) {
    std::string rMaskExpr = actionArgs.GetMaskNext();
    if (rMaskExpr.empty())
      rMaskExpr = Mask_.MaskExpression();
    refMask.SetMaskString( rMaskExpr );
    if (refFrm.Parm().SetupIntegerMask( refMask, refFrm.Coord() ))
      return Action::ERR;
    // Get center of mask in reference
    if (useMass_)
      refCenter_ = refFrm.Coord().VCenterOfMass( refMask );
    else
      refCenter_ = refFrm.Coord().VGeometricCenter( refMask );
    centerMode_ = Frame::POINT; 
  }

  mprintf("    CENTER: Centering coordinates using");
  if (useMass_)
    mprintf(" center of mass");
  else
    mprintf(" geometric center");
  mprintf(" of atoms in mask (%s) to\n", Mask_.MaskString());
  if (centerMode_ == Frame::POINT)
    mprintf("\tcenter of mask (%s) in reference '%s'.\n", refMask.MaskString(),
            refFrm.FrameName().base());
  else if (centerMode_ == Frame::ORIGIN)
    mprintf("\tcoordinate origin.\n");
  else
    mprintf("\tbox center.\n");

  return Action::OK;
}

// Action_Center::Setup()
Action::RetType Action_Center::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask(Mask_) ) return Action::ERR;
  Mask_.MaskInfo();
  if (Mask_.None()) {
    mprintf("Warning: Mask contains 0 atoms.\n");
    return Action::ERR;
  }

  if (centerMode_ == Frame::BOXCTR && currentParm->BoxType()==Box::NOBOX) {
    mprintf("Warning: Box center specified but no box information.\n");
    //fprintf(stdout,"                            Centering on origin.\n");
    return Action::ERR;
  }

  return Action::OK;  
}

// Action_Center::DoAction()
/** Center coordinates in frame according to specified mode. */
Action::RetType Action_Center::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  currentFrame->Center(Mask_, centerMode_, refCenter_, useMass_);

  return Action::OK;
}
