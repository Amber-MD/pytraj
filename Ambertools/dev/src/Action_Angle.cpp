// Action_Angle 
#include "Action_Angle.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_Angle::Action_Angle() :
  ang_(0),
  useMass_(false)
{ } 

void Action_Angle::Help() {
  mprintf("\t[<name>] <mask1> <mask2> <mask3> [out <filename>] [mass]\n"
          "  Calculate the angle between atoms in masks 1-3.\n");
}

// Action_Angle::init()
Action::RetType Action_Angle::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  useMass_ = actionArgs.hasKey("mass");

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  std::string mask2 = actionArgs.GetMaskNext();
  std::string mask3 = actionArgs.GetMaskNext();
  if (mask1.empty() || mask2.empty() || mask3.empty()) {
    mprinterr("Error: angle: Requires 3 masks\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);
  Mask2_.SetMaskString(mask2);
  Mask3_.SetMaskString(mask3);

  // Dataset to store angles
  ang_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"Ang");
  if (ang_==0) return Action::ERR;
  ang_->SetScalar( DataSet::M_ANGLE );
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( ang_ );

  mprintf("    ANGLE: [%s]-[%s]-[%s]\n",Mask1_.MaskString(), Mask2_.MaskString(), 
          Mask3_.MaskString());
  if (useMass_)
    mprintf("              Using center of mass of atoms in masks.\n");

  return Action::OK;
}

// Action_Angle::setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_Angle::Setup(Topology* currentParm, Topology** parmAddress) {

  if (currentParm->SetupIntegerMask(Mask1_)) return Action::ERR;
  if (currentParm->SetupIntegerMask(Mask2_)) return Action::ERR;
  if (currentParm->SetupIntegerMask(Mask3_)) return Action::ERR;
  mprintf("\t");
  Mask1_.BriefMaskInfo();
  Mask2_.BriefMaskInfo();
  Mask3_.BriefMaskInfo();
  mprintf("\n");
  if (Mask1_.None() || Mask2_.None() || Mask3_.None()) {
    mprintf("Warning: angle: One or more masks contain 0 atoms.\n");
    return Action::ERR;
  }

  return Action::OK;  
}

// Action_Angle::action()
Action::RetType Action_Angle::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Vec3 a1, a2, a3;
  if (useMass_) {
    a1 = currentFrame->VCenterOfMass( Mask1_ );
    a2 = currentFrame->VCenterOfMass( Mask2_ );
    a3 = currentFrame->VCenterOfMass( Mask3_ );
  } else {
    a1 = currentFrame->VGeometricCenter( Mask1_ );
    a2 = currentFrame->VGeometricCenter( Mask2_ );
    a3 = currentFrame->VGeometricCenter( Mask3_ );
  }
  double aval = CalcAngle( a1.Dptr(), a2.Dptr(), a3.Dptr() );

  aval *= Constants::RADDEG;

  ang_->Add(frameNum, &aval);

  return Action::OK;
}

