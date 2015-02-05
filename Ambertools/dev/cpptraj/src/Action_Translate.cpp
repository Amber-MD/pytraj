#include "Action_Translate.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Translate::Action_Translate() { }

void Action_Translate::Help() {
  mprintf("\t[<mask>] [x <dx>] [y <dy>] [z <dz>]\n"
          "  Translate atoms in <mask>\n");
}

Action::RetType Action_Translate::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  double x = actionArgs.getKeyDouble("x",0.0);
  double y = actionArgs.getKeyDouble("y",0.0);
  double z = actionArgs.getKeyDouble("z",0.0);
  Trans_.SetVec(x, y, z);
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    TRANSLATE: Translating atoms in mask %s\n", mask_.MaskString());
  mprintf("\t%f Ang. in X, %f Ang. in Y, %f Ang. in Z\n",
          Trans_[0], Trans_[1], Trans_[2]);
  return Action::OK;
};

Action::RetType Action_Translate::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: translate: No atoms selected.\n");
    return Action::ERR;
  }
  return Action::OK;
}

Action::RetType Action_Translate::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    currentFrame->Translate(Trans_, *atom);
  return Action::OK;
}
