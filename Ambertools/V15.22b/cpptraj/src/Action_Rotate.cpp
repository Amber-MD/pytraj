#include "Action_Rotate.h"
#include "CpptrajStdio.h"
#include "Constants.h"

// CONSTRUCTOR
Action_Rotate::Action_Rotate() { }

void Action_Rotate::Help() {
  mprintf("\t[<mask>] [x <xdeg>] [y <ydeg>] [z <zdeg>]\n"
          "  Rotate atoms in <mask> around x, y, and/or z axes.\n");
}

Action::RetType Action_Rotate::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  double xrot = actionArgs.getKeyDouble("x",0.0);
  double yrot = actionArgs.getKeyDouble("y",0.0);
  double zrot = actionArgs.getKeyDouble("z",0.0);
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  // Calc rotation matrix
  RotMatrix_.CalcRotationMatrix( xrot * Constants::DEGRAD, yrot * Constants::DEGRAD, 
                                 zrot * Constants::DEGRAD );

  mprintf("    ROTATE: Rotating atoms in mask %s\n", mask_.MaskString());
  mprintf("\t%f degrees around X, %f degrees around Y, %f degrees around Z\n",
          xrot, yrot, zrot);
  return Action::OK;
};

Action::RetType Action_Rotate::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: rotate: No atoms selected.\n");
    return Action::ERR;
  }
  return Action::OK;
}

Action::RetType Action_Rotate::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  currentFrame->Rotate(RotMatrix_, mask_);
  return Action::OK;
}

