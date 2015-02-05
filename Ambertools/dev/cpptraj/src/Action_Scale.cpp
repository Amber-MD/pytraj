#include "Action_Scale.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Scale::Action_Scale() :
  sx_(1),
  sy_(1),
  sz_(1)
{}

void Action_Scale::Help() {
  mprintf("\t[x <sx>] [y <sy>] [z <sz>] [<mask>]\n"
          "\tScale the position of atoms in <mask>\n");
}

// Action_Scale::init()
Action::RetType Action_Scale::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  sx_ = actionArgs.getKeyDouble("x", 1);
  sy_ = actionArgs.getKeyDouble("y", 1);
  sz_ = actionArgs.getKeyDouble("z", 1);
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    SCALE coordinates: X by %.3f, Y by %.3f, Z by %.3f\n", sx_, sy_, sz_);
  mprintf("                       Mask is [%s]\n", mask_.MaskString());

  return Action::OK;
}

// Action_Scale::setup()
Action::RetType Action_Scale::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return Action::ERR;
  if ( mask_.None() ) {
    mprintf("Warning: scale: No atoms selected.\n");
    return Action::ERR;
  }
  return Action::OK;
}

// Action_Scale::action()
Action::RetType Action_Scale::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  currentFrame->Scale(mask_, sx_, sy_, sz_);
  return Action::OK;
}
