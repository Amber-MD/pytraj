#include "Action_Principal.h"
#include "CpptrajStdio.h"
#include "Matrix_3x3.h"

// CONSTRUCTOR
Action_Principal::Action_Principal() :
  doRotation_(false),
  useMass_(false),
  debug_(0)
{ }

void Action_Principal::Help() {
  mprintf("\t[<mask>] [dorotation] [mass] [out <filename>]\n"
          "  Calculate principal axes of atoms in <mask>. Align the system along\n"
          "  principal axes if 'dorotation' specified.\n");
}

// Action_Principal::init()
Action::RetType Action_Principal::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Keywords
  doRotation_ = actionArgs.hasKey("dorotation");
  useMass_ = actionArgs.hasKey("mass");
  std::string filename = actionArgs.GetStringKey("out");
  if (!doRotation_ && filename.empty()) {
    mprinterr("Error: principal: At least one of 'dorotation' or 'out <filename>' must be specified.\n");
    return Action::ERR;
  }
  // Masks
  mask_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    PRINCIPAL:");
  if (!filename.empty()) {
    mprintf(" output eigenvectors/eigenvalues to %s,", filename.c_str());
    if (outfile_.OpenEnsembleWrite(filename, DSL->EnsembleNum())) return Action::ERR;
  }
  if (doRotation_)
    mprintf(" with rotation by");
  else
    mprintf(" without rotation by");
  if (useMass_)
    mprintf(" center of mass");
  else
    mprintf(" center of geometry");
  mprintf(", atoms selected by [%s]\n", mask_.MaskString());

  return Action::OK;
}

// Action_Principal::setup()
Action::RetType Action_Principal::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask(mask_)) return Action::ERR;

  if (mask_.None()) {
    mprintf("Warning: No atoms selected for %s [%s].\n",currentParm->c_str(), mask_.MaskString());
    return Action::ERR;
  }

  mprintf("\tSelected %i atoms.\n", mask_.Nselected());
  return Action::OK;
}

// Action_Principal::action()
Action::RetType Action_Principal::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 Inertia;
  Vec3 Eval;

  currentFrame->CalculateInertia( mask_, Inertia );

  // NOTE: Diagonalize_Sort_Chirality places sorted eigenvectors in rows.
  Inertia.Diagonalize_Sort_Chirality( Eval, debug_ );
  if (outfile_.IsOpen()) {
    int fn = frameNum+1; 
    outfile_.Printf("%i EIGENVALUES: %f %f %f\n%i EIGENVECTOR 0: %f %f %f\n%i EIGENVECTOR 1: %f %f %f\n%i EIGENVECTOR 2: %f %f %f\n", 
      fn, Eval[0], Eval[1], Eval[2],
      fn, Inertia[0], Inertia[1], Inertia[2],
      fn, Inertia[3], Inertia[4], Inertia[5],
      fn, Inertia[6], Inertia[7], Inertia[8]);
    //Eval.Print("PRINCIPAL EIGENVALUES");
    //Inertia.Print("PRINCIPAL EIGENVECTORS (Rows)");
  }
  
  // Rotate - since Evec is already transposed (eigenvectors
  // are returned in rows) just do plain rotation to affect an
  // inverse rotation.
  if (doRotation_)
    currentFrame->Rotate( Inertia );

  return Action::OK;
}
