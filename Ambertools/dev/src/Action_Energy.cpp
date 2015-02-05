#include "Action_Energy.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Energy::Action_Energy() :
  currentParm_(0)
{}


void Action_Energy::Help() {
  mprintf("\t[<name>] [<mask1>] [out <filename>]\n"
          "\t[bond] [angle] [dihedral] [nb14] [nonbond]\n"
          "  Calculate energy for atoms in mask.\n");
}

/// DataSet aspects
static const char* Estring[] = {"bond", "angle", "dih", "vdw14", "elec14", "vdw", "elec", "total"};

/// Calculation types
static const char* Cstring[] = {"Bond", "Angle", "Torsion", "1-4_Nonbond", "Nonbond" };

int Action_Energy::AddSet(Etype typeIn, DataSetList* DSL, DataFile* outfile,
                          std::string const& setname)
{
  Energy_[typeIn] = DSL->AddSetAspect(DataSet::DOUBLE, setname, Estring[typeIn]);
  if (Energy_[typeIn] == 0) return 1;
  if (outfile != 0) outfile->AddSet( Energy_[typeIn] );
  return 0;
}

// Action_Energy::Init()
Action::RetType Action_Energy::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ENE_.SetDebug( debugIn );
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Which terms will be calculated?
  Ecalcs_.clear();
  if (actionArgs.hasKey("bond"))     Ecalcs_.push_back(BND);
  if (actionArgs.hasKey("angle"))    Ecalcs_.push_back(ANG);
  if (actionArgs.hasKey("dihedral")) Ecalcs_.push_back(DIH);
  if (actionArgs.hasKey("nb14"))     Ecalcs_.push_back(N14);
  if (actionArgs.hasKey("nonbond"))  Ecalcs_.push_back(NBD);
  // If nothing is selected, select all.
  if (Ecalcs_.empty()) {
    for (int c = 0; c <= (int)NBD; c++)
      Ecalcs_.push_back( (CalcType)c );
  }

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // DataSet
  std::string setname = actionArgs.GetStringNext();
  if (setname.empty())
    setname = DSL->GenerateDefaultName("ENE");
  Energy_.clear();
  Energy_.resize( (int)TOTAL + 1, 0 );
  for (calc_it calc = Ecalcs_.begin(); calc != Ecalcs_.end(); ++calc)
  {
    switch (*calc) {
      case BND: if (AddSet(BOND, DSL, outfile, setname)) return Action::ERR; break;
      case ANG: if (AddSet(ANGLE, DSL, outfile, setname)) return Action::ERR; break;
      case DIH: if (AddSet(DIHEDRAL, DSL, outfile, setname)) return Action::ERR; break;
      case N14:
        if (AddSet(V14, DSL, outfile, setname)) return Action::ERR;
        if (AddSet(Q14, DSL, outfile, setname)) return Action::ERR;
        break;
      case NBD:
        if (AddSet(VDW, DSL, outfile, setname)) return Action::ERR;
        if (AddSet(ELEC, DSL, outfile, setname)) return Action::ERR;
        break;
    }
  }
//  if (Ecalcs_.size() > 1) {
    if (AddSet(TOTAL, DSL, outfile, setname)) return Action::ERR;
//  }
      
  mprintf("    ENERGY: Calculating energy for atoms in mask '%s'\n", Mask1_.MaskString());
  mprintf("\tCalculating terms:");
  for (calc_it calc = Ecalcs_.begin(); calc != Ecalcs_.end(); ++calc)
    mprintf(" %s", Cstring[*calc]);
  mprintf("\n");

  return Action::OK;
}

// Action_Energy::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_Energy::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupCharMask(Mask1_)) return Action::ERR;
  if (Mask1_.None()) {
    mprinterr("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
    return Action::ERR;
  }
  Mask1_.MaskInfo();
  Imask_ = Mask1_;
  Imask_.ConvertToIntMask();
  currentParm_ = currentParm;
  return Action::OK;
}

// Action_Energy::DoAction()
Action::RetType Action_Energy::DoAction(int frameNum, Frame* currentFrame,
                                            Frame** frameAddress)
{
  double Etot = 0.0, ene, ene2;
  for (calc_it calc = Ecalcs_.begin(); calc != Ecalcs_.end(); ++calc)
  {
    switch (*calc) {
      case BND:
        ene = ENE_.E_bond(*currentFrame, *currentParm_, Mask1_);
        Energy_[BOND]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case ANG:
        ene = ENE_.E_angle(*currentFrame, *currentParm_, Mask1_);
        Energy_[ANGLE]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case DIH:
        ene = ENE_.E_torsion(*currentFrame, *currentParm_, Mask1_);
        Energy_[DIHEDRAL]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case N14:
        ene = ENE_.E_14_Nonbond(*currentFrame, *currentParm_, Mask1_, ene2);
        Energy_[V14]->Add(frameNum, &ene);
        Energy_[Q14]->Add(frameNum, &ene2);
        Etot += (ene + ene2);
        break;
      case NBD:
        ene = ENE_.E_Nonbond(*currentFrame, *currentParm_, Imask_, ene2);
        Energy_[VDW]->Add(frameNum, &ene);
        Energy_[ELEC]->Add(frameNum, &ene2);
        Etot += (ene + ene2);
        break;
    }
  }

  Energy_[TOTAL]->Add(frameNum, &Etot);

  return Action::OK;
}

void Action_Energy::Print() {
  mprintf("Timing for energy: '%s' ('%s')\n", Energy_[TOTAL]->Legend().c_str(),
           Mask1_.MaskString());
  ENE_.PrintTiming();
} 
