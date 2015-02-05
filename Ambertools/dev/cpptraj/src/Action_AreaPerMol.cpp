#include "Action_AreaPerMol.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_AreaPerMol::Action_AreaPerMol() :
  area_per_mol_(0),
  Nmols_(-1.0),
  Nlayers_(1.0),
  areaType_(XY)
{}

void Action_AreaPerMol::Help() {
  mprintf("\t[<name>] {[<mask1>] [nlayers <#>] | nmols <#>} [out <filename>] [{xy | xz | yz}]\n"
          "  Calculate the specified area per molecule for molecules in <mask1>.\n");
}

static const char* APMSTRING[] = {"XY", "XZ", "YZ"};

// Action_AreaPerMol::Init()
Action::RetType Action_AreaPerMol::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  if (actionArgs.hasKey("xy")) areaType_ = XY;
  else if (actionArgs.hasKey("xz")) areaType_ = XZ;
  else if (actionArgs.hasKey("yz")) areaType_ = YZ;
  else areaType_ = XY;

  Nmols_ = (double)actionArgs.getKeyInt("nmols", -1);

  // Get Masks
  if (Nmols_ < 0.0) {
    Nlayers_ = (double)actionArgs.getKeyInt("nlayers", 1);
    if (Nlayers_ < 1.0) {
      mprinterr("Error: Number of layers must be > 0\n");
      return Action::ERR;
    }
    Mask1_.SetMaskString( actionArgs.GetMaskNext() );
  }

  // DataSet
  area_per_mol_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"APM");
  if (area_per_mol_==0) return Action::ERR;
  // Add DataSet to DataFileList
  if (outfile != 0) outfile->AddSet( area_per_mol_ );

  mprintf("    AREAPERMOL: Calculating %s area per molecule", APMSTRING[areaType_]);
  if (Mask1_.MaskStringSet())
    mprintf(" using mask '%s', %.0f layers.\n", Mask1_.MaskString(), Nlayers_);
  else
    mprintf(" for %.0f mols\n", Nmols_);

  return Action::OK;
}

// Action_AreaPerMol::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_AreaPerMol::Setup(Topology* currentParm, Topology** parmAddress) {
  // Needs box info
  if (currentParm->BoxType() == Box::NOBOX) {
    mprintf("Warning: No box information for '%s', cannot calculate area.\n",
            currentParm->c_str());
    return Action::ERR;
  }
  // Probably will not work for non-orthorhombic cells
  if (currentParm->BoxType() != Box::ORTHO)
    mprintf("Warning: Box is not orthorhombic, calculated area may not be correct.\n");
  // Determine how many molecules are selected
  if (Mask1_.MaskStringSet()) {
    if (currentParm->SetupCharMask(Mask1_)) return Action::ERR;
    if (Mask1_.None()) {
      mprinterr("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
      return Action::ERR;
    }
    Nmols_ = 0.0;
    for (Topology::mol_iterator mol = currentParm->MolStart();
                                mol != currentParm->MolEnd(); ++mol)
    {
      if (Mask1_.AtomsInCharMask(mol->BeginAtom(), mol->EndAtom()))
        Nmols_ += 1.0;
    }
    mprintf("\tMask '%s' selects %.0f molecules.\n", Mask1_.MaskString(), Nmols_);
    if (Nmols_ < 1.0) return Action::ERR;
    Nmols_ /= Nlayers_;
    mprintf("\tArea per %.0f molecules (%0.f layers) will be determined.\n", Nmols_, Nlayers_);
  }
  return Action::OK;
}

// Action_AreaPerMol::DoAction()
Action::RetType Action_AreaPerMol::DoAction(int frameNum, Frame* currentFrame,
                                            Frame** frameAddress)
{
  double area;
  if (areaType_ == XY)
    area = currentFrame->BoxCrd().BoxX() * currentFrame->BoxCrd().BoxY();
  else if (areaType_ == XZ) 
    area = currentFrame->BoxCrd().BoxX() * currentFrame->BoxCrd().BoxZ();
  else // if areaType_ == YZ
    area = currentFrame->BoxCrd().BoxY() * currentFrame->BoxCrd().BoxZ();

  area = area / Nmols_;

  area_per_mol_->Add(frameNum, &area);

  return Action::OK;
}
