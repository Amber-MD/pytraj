#include "Action_MultiDihedral.h"
#include "DataSet_double.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "TorsionRoutines.h"
#include "StringRoutines.h" // convertToInteger

Action_MultiDihedral::Action_MultiDihedral() :
  debug_(0),
  range360_(false),
  outfile_(0),
  masterDSL_(0)
{}

void Action_MultiDihedral::Help() {
  mprintf("\t[<name>] <dihedral types> [resrange <range>] [out <filename>] [range360]\n");
  mprintf("\t[dihtype <name>:<a0>:<a1>:<a2>:<a3>[:<offset>] ...]\n");
  DihedralSearch::OffsetHelp();
  //mprintf("\t[range360]\n");
  mprintf("\t<dihedral types> = ");
  DihedralSearch::ListKnownTypes();
  mprintf("  Calculate specified dihedral angle types for residues in given <range>.\n");
}

Action::RetType Action_MultiDihedral::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  outfile_ = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  range360_ = actionArgs.hasKey("range360");
  std::string resrange_arg = actionArgs.GetStringKey("resrange");
  if (!resrange_arg.empty())
    if (resRange_.SetRange( resrange_arg )) return Action::ERR;
  // Search for known dihedral keywords
  dihSearch_.SearchForArgs(actionArgs);
  // Get custom dihedral arguments: dihtype <name>:<a0>:<a1>:<a2>:<a3>[:<offset>]
  std::string dihtype_arg = actionArgs.GetStringKey("dihtype");
  while (!dihtype_arg.empty()) {
    ArgList dihtype(dihtype_arg, ":");
    if (dihtype.Nargs() < 5) {
      mprinterr("Error: Malformed dihtype arg.\n");
      return Action::ERR;
    }
    int offset = 0;
    if (dihtype.Nargs() == 6) offset = convertToInteger(dihtype[5]);
    dihSearch_.SearchForNewType(offset,dihtype[1],dihtype[2],dihtype[3],dihtype[4], dihtype[0]);
    dihtype_arg = actionArgs.GetStringKey("dihtype");
  }
  // If no dihedral types yet selected, this will select all.
  dihSearch_.SearchForAll();

  // Setup DataSet(s) name
  dsetname_ = actionArgs.GetStringNext();

  mprintf("    MULTIDIHEDRAL: Calculating");
  dihSearch_.PrintTypes();
  if (!resRange_.Empty())
    mprintf(" dihedrals for residues in range %s\n", resRange_.RangeArg());
  else
    mprintf(" dihedrals for all solute residues.\n");
  if (!dsetname_.empty())
    mprintf("\tDataSet name: %s\n", dsetname_.c_str());
  if (outfile_ != 0) mprintf("\tOutput to %s\n", outfile_->DataFilename().base());
  if (range360_) 
    mprintf("\tOutput range is 0 to 360 degrees.\n");
  else
    mprintf("\tOutput range is -180 to 180 degrees.\n");
  DSL->SetDataSetsPending(true);
  masterDSL_ = DSL;
  return Action::OK;
}

// Action_MultiDihedral::Setup();
Action::RetType Action_MultiDihedral::Setup(Topology* currentParm, Topology** parmAddress) {
  Range actualRange;
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  if (resRange_.Empty())
    actualRange = currentParm->SoluteResidues();
  else {
    // If user range specified, create new range shifted by -1 since internal
    // resnums start from 0.
    actualRange = resRange_;
    actualRange.ShiftBy(-1);
  }
  // Exit if no residues specified
  if (actualRange.Empty()) {
    mprinterr("Error: No residues specified for %s\n",currentParm->c_str());
    return Action::ERR;
  }
  // Search for specified dihedrals in each residue in the range
  if (dihSearch_.FindDihedrals(*currentParm, actualRange))
    return Action::ERR;
  mprintf("\tResRange=[%s]", resRange_.RangeArg());
  dihSearch_.PrintTypes();
  mprintf(", %i dihedrals.\n", dihSearch_.Ndihedrals());

  // Print selected dihedrals, set up DataSets
  data_.clear();
  if (dsetname_.empty())
    dsetname_ = masterDSL_->GenerateDefaultName("MDIH");
  for (DihedralSearch::mask_it dih = dihSearch_.begin();
                               dih != dihSearch_.end(); ++dih)
  {
    int resNum = dih->ResNum() + 1;
    // See if Dataset already present
    DataSet* ds = masterDSL_->GetSet(dsetname_, resNum, dih->Name());
    if (ds == 0) {
      // Create new DataSet
      ds = masterDSL_->AddSetIdxAspect( DataSet::DOUBLE, dsetname_, resNum, dih->Name());
      if (ds == 0) return Action::ERR;
      // FIXME: Dihedral types in DihedralSearch.h and DataSet.h should be
      //        consolidated.
      DataSet::scalarType dstype = DataSet::UNDEFINED;
      switch (dih->Type()) {
        case DihedralSearch::PHI: dstype = DataSet::PHI; break;
        case DihedralSearch::PSI: dstype = DataSet::PSI; break;
        case DihedralSearch::CHIP: dstype = DataSet::PCHI; break;
        case DihedralSearch::ALPHA: dstype = DataSet::ALPHA; break;
        case DihedralSearch::BETA: dstype = DataSet::BETA; break;
        case DihedralSearch::GAMMA: dstype = DataSet::GAMMA; break;
        case DihedralSearch::DELTA: dstype = DataSet::DELTA; break;
        case DihedralSearch::EPSILON: dstype = DataSet::EPSILON; break;
        case DihedralSearch::ZETA: dstype = DataSet::ZETA; break;
        case DihedralSearch::CHIN: dstype = DataSet::CHI; break;
        case DihedralSearch::OMEGA: dstype = DataSet::OMEGA; break;
        default: dstype = DataSet::UNDEFINED;
      }
      ds->SetScalar( DataSet::M_TORSION, dstype );
      // Add to outfile
      if (outfile_ != 0)
        outfile_->AddSet( ds );
    }
    data_.push_back( ds ); 
    if (debug_ > 0) {
      mprintf("\tDIH [%s]:", ds->Legend().c_str());
      mprintf(" :%i@%i",   (*currentParm)[dih->A0()].ResNum()+1, dih->A0() + 1);
      mprintf(" :%i@%i",   (*currentParm)[dih->A1()].ResNum()+1, dih->A1() + 1);
      mprintf(" :%i@%i",   (*currentParm)[dih->A2()].ResNum()+1, dih->A2() + 1);
      mprintf(" :%i@%i\n", (*currentParm)[dih->A3()].ResNum()+1, dih->A3() + 1);
    }
  }
  return Action::OK;
}

// Action_MultiDihedral::DoAction()
Action::RetType Action_MultiDihedral::DoAction(int frameNum, Frame* currentFrame, 
                                               Frame** frameAddress)
{
  std::vector<DataSet*>::const_iterator ds = data_.begin();
  for (DihedralSearch::mask_it dih = dihSearch_.begin();
                               dih != dihSearch_.end(); ++dih, ++ds)
  {
    double torsion = Torsion( currentFrame->XYZ(dih->A0()),
                              currentFrame->XYZ(dih->A1()),
                              currentFrame->XYZ(dih->A2()),
                              currentFrame->XYZ(dih->A3()) );
    torsion *= Constants::RADDEG;
    (*ds)->Add(frameNum, &torsion);
  }
  return Action::OK;
}

void Action_MultiDihedral::Print() {
  if (range360_) {
    for (std::vector<DataSet*>::const_iterator ds = data_.begin();
                                               ds != data_.end(); ++ds)
      ((DataSet_double*)*ds)->ShiftTorsions(0.0, 0.0);
  }
}
