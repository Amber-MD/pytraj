#include "Action_Temperature.h"
#include "CpptrajStdio.h"

Action_Temperature::Action_Temperature() :
  Tdata_(0),
  getTempFromFrame_(false),
  shakeType_(OFF),
  degrees_of_freedom_(0)
{}

void Action_Temperature::Help() {
  mprintf("\t[<name>] {frame | [<mask>] [ntc <#>]} [out <filename>]\n"
          "  Calculate temperature in frame based on velocity information.\n"
          "  If 'frame' is specified just use frame temperature (read in from\n"
          "  e.g. REMD trajectory)\n");
}

static const char* ShakeString[] = {
  "off", "bonds to H", "all bonds"
};

// Action_Temperature::Init()
Action::RetType Action_Temperature::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Keywords
  if (actionArgs.hasKey("frame")) {
    getTempFromFrame_ = true;
    shakeType_ = OFF;
    degrees_of_freedom_ = 0;
  } else {
    getTempFromFrame_ = false;
    int ntc = actionArgs.getKeyInt("ntc",-1);
    if (ntc != -1) {
      if (ntc < 1 || ntc > 3) {
        mprinterr("Error: temperature: ntc must be 1, 2, or 3\n");
        return Action::ERR;
      }
      shakeType_ = (ShakeType)(ntc - 1);
    } else
      shakeType_ = OFF;
  }
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Masks
  if (!getTempFromFrame_)
    Mask_.SetMaskString( actionArgs.GetMaskNext() );
  // DataSet 
  Tdata_ =  DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "Tdata");
  if (Tdata_ == 0) return Action::ERR;
  if (outfile != 0) outfile->AddSet( Tdata_ );
  
  if (getTempFromFrame_) {
    mprintf("    TEMPERATURE: Frame temperatures will be saved in data set %s\n",
             Tdata_->Legend().c_str());
  } else {
    mprintf("    TEMPERATURE: Calculate temperature for atoms in mask [%s]\n", Mask_.MaskString());
    mprintf("\tUsing SHAKE (ntc) value of [%s]\n", ShakeString[shakeType_]);
  }
  return Action::OK;
}

// Action_Temperature::Setup()
Action::RetType Action_Temperature::Setup(Topology* currentParm, Topology** parmAddress)
{
  if (!getTempFromFrame_) {
    // Masks
    if (currentParm->SetupIntegerMask( Mask_ )) return Action::ERR;
    Mask_.MaskInfo();
    if (Mask_.None()) {
      mprintf("Warning: temperature: No atoms selected in [%s]\n", Mask_.MaskString());
      return Action::ERR;
    }
    // Calculate degrees of freedom
    // If SHAKE is on, add up all bonds which cannot move because of SHAKE.
    // NOTE: For now, dont distinguish between solute/solvent.
    int constrained_bonds_to_h = 0;
    int constrained_heavy_bonds = 0;
    if (shakeType_ >= BONDS_TO_H) {
      constrained_bonds_to_h = (int)currentParm->BondsH().size();
      mprintf("\t%i bonds to hydrogen constrained.\n", constrained_bonds_to_h);
      if (shakeType_ >= ALL_BONDS) {
        constrained_heavy_bonds = (int)currentParm->Bonds().size();
        mprintf("\t%i bonds to heavy atoms constrained.\n", constrained_heavy_bonds);
      }
    }
    // Just estimate for now, 3N - 6
    degrees_of_freedom_ = (3 * Mask_.Nselected()) - constrained_bonds_to_h
                          - constrained_heavy_bonds - 6;
    mprintf("\t# of degrees of freedom = %i\n", degrees_of_freedom_);
  }
  return Action::OK;
}

// Action_Temperature::DoAction()
Action::RetType Action_Temperature::DoAction(int frameNum, Frame* currentFrame, 
                                             Frame** frameAddress) 
{
  double tdata;
  if (getTempFromFrame_)
    tdata = currentFrame->Temperature();
  else
    tdata = currentFrame->CalcTemperature(Mask_, degrees_of_freedom_);
  Tdata_->Add(frameNum, &tdata);
  return Action::OK;
}
