// Action_Pucker
#include "Action_Pucker.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "TorsionRoutines.h"
#include "DataSet_double.h"

// CONSTRUCTOR
Action_Pucker::Action_Pucker() :
  pucker_(0),
  amplitude_(0),
  theta_(0),
  puckerMethod_(ALTONA),
  useMass_(true),
  range360_(false),
  offset_(0.0)
{ } 

void Action_Pucker::Help() {
  mprintf("\t[<name>] <mask1> <mask2> <mask3> <mask4> <mask5> [<mask6>] [geom]\n"
          "\t[out <filename>] [altona | cremer] [amplitude] [theta]\n"
          "\t[range360] [offset <offset>]\n"
          "\tCalculate pucker of atoms in masks 1-5 (or 6, 'cremer' only).\n");
}

// Action_Pucker::Init()
Action::RetType Action_Pucker::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  if      (actionArgs.hasKey("altona")) puckerMethod_=ALTONA;
  else if (actionArgs.hasKey("cremer")) puckerMethod_=CREMER;
  bool calc_amp = actionArgs.hasKey("amplitude");
  bool calc_theta = actionArgs.hasKey("theta");
  offset_ = actionArgs.getKeyDouble("offset",0.0);
  range360_ = actionArgs.hasKey("range360");
  useMass_ = !actionArgs.hasKey("geom");
  DataSet::scalarType stype = DataSet::UNDEFINED;
  std::string stypename = actionArgs.GetStringKey("type");
  if ( stypename == "pucker" ) stype = DataSet::PUCKER;

  // Get Masks
  Masks_.clear();
  std::string mask_expression = actionArgs.GetMaskNext();
  while (!mask_expression.empty()) {
    Masks_.push_back( AtomMask( mask_expression ) );
    mask_expression = actionArgs.GetMaskNext();
  }
  if (Masks_.size() < 5 || Masks_.size() > 6) {
    mprinterr("Error: Pucker can only be calculated for 5 or 6 masks, %zu specified.\n",
              Masks_.size());
    return Action::ERR;
  }
  if (Masks_.size() > 5 && puckerMethod_ != CREMER) {
    mprinterr("Error: Pucker with %zu masks only supported with 'cremer'\n");
    return Action::ERR;
  }
  // Set up array to hold coordinate vectors.
  AX_.resize( Masks_.size() );

  // Setup dataset
  pucker_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"Pucker");
  if (pucker_ == 0) return Action::ERR;
  pucker_->SetScalar( DataSet::M_PUCKER, stype );
  amplitude_ = 0;
  theta_ = 0;
  if (calc_amp)
    amplitude_ = DSL->AddSetAspect(DataSet::DOUBLE, pucker_->Name(), "Amp");
  if (calc_theta) {
    if ( Masks_.size() < 6 )
      mprintf("Warning: 'theta' calc. not supported for < 6 masks.\n");
    else
      theta_ = DSL->AddSetAspect(DataSet::DOUBLE, pucker_->Name(), "Theta");
  }
  // Add dataset to datafile list
  if (outfile != 0) {
    outfile->AddSet( pucker_ );
    if (amplitude_ != 0) outfile->AddSet( amplitude_ );
    if (theta_ != 0) outfile->AddSet( theta_ );
  }

  mprintf("    PUCKER: ");
  for (std::vector<AtomMask>::const_iterator MX = Masks_.begin();
                                             MX != Masks_.end(); ++MX)
  {
    if (MX != Masks_.begin()) mprintf("-");
    mprintf("[%s]", MX->MaskString());
  }
  mprintf("\n");
  if (puckerMethod_==ALTONA) 
    mprintf("            Using Altona & Sundaralingam method.\n");
  else if (puckerMethod_==CREMER)
    mprintf("            Using Cremer & Pople method.\n");
  if (outfile != 0) 
    mprintf("            Data will be written to %s\n", outfile->DataFilename().base());
  if (amplitude_!=0)
    mprintf("            Amplitudes will be stored.\n");
  if (theta_!=0)
    mprintf("            Thetas will be stored.\n");
  if (offset_!=0)
    mprintf("            Offset: %f will be added to values.\n");
  if (range360_)
    mprintf("              Output range is 0 to 360 degrees.\n");
  else
    mprintf("              Output range is -180 to 180 degrees.\n");

  return Action::OK;
}

// Action_Pucker::Setup()
Action::RetType Action_Pucker::Setup(Topology* currentParm, Topology** parmAddress) {
  mprintf("\t");
  for (std::vector<AtomMask>::iterator MX = Masks_.begin();
                                       MX != Masks_.end(); ++MX)
  {
    if ( currentParm->SetupIntegerMask( *MX ) ) return Action::ERR;
    MX->BriefMaskInfo();
    if (MX->None()) {
      mprintf("\nWarning: Mask '%s' selects no atoms for topology '%s'\n",
              MX->MaskString(), currentParm->c_str());
      return Action::ERR;
    }
  }
  mprintf("\n");

  return Action::OK;  
}

// Action_Pucker::DoAction()
Action::RetType Action_Pucker::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double pval, aval, tval;
  std::vector<Vec3>::iterator ax = AX_.begin(); 

  if (useMass_) {
    for (std::vector<AtomMask>::const_iterator MX = Masks_.begin();
                                               MX != Masks_.end(); ++MX, ++ax)
      *ax = currentFrame->VCenterOfMass( *MX );
  } else {
     for (std::vector<AtomMask>::const_iterator MX = Masks_.begin();
                                               MX != Masks_.end(); ++MX, ++ax)
      *ax = currentFrame->VGeometricCenter( *MX );
  }

  switch (puckerMethod_) {
    case ALTONA: 
      pval = Pucker_AS( AX_[0].Dptr(), AX_[1].Dptr(), AX_[2].Dptr(), 
                        AX_[3].Dptr(), AX_[4].Dptr(), aval );
      break;
    case CREMER:
      pval = Pucker_CP( AX_[0].Dptr(), AX_[1].Dptr(), AX_[2].Dptr(), 
                        AX_[3].Dptr(), AX_[4].Dptr(), AX_[5].Dptr(), 
                        AX_.size(), aval, tval );
      break;
  }
  if ( amplitude_ != 0 )
    amplitude_->Add(frameNum, &aval);
  if ( theta_ != 0 ) {
    tval *= Constants::RADDEG;
    theta_->Add(frameNum, &tval);
  }
  pval *= Constants::RADDEG;
  pucker_->Add(frameNum, &pval);

  return Action::OK;
} 

void Action_Pucker::Print() {
  double puckermin;
  if (range360_)
    puckermin = 0.0;
  else
    puckermin = -180.0;
  ((DataSet_double*)pucker_)->ShiftTorsions(puckermin, offset_);
}
