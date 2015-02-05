#include "Action_CheckChirality.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_CheckChirality::Action_CheckChirality() :
  ensembleNum_(-1)
  //outfile_(0), masterDSL_(0)
{}

void Action_CheckChirality::Help() {
  mprintf("\t[<name>] [<mask1>] [out <filename>]\n"
          "  Check the chirality of AA residues in <mask1>.\n");
}

// Action_CheckChirality::Init()
Action::RetType Action_CheckChirality::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  // Get keywords
  //outfile_ = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  outfilename_ = actionArgs.GetStringKey("out");

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // DataSet name
  //setname_ = actionArgs.GetStringNext();

  mprintf("    CHECKCHIRALITY: Check chirality for AA residues in mask '%s'\n",
          Mask1_.MaskString());
  if (!outfilename_.empty())
    mprintf("\tOutput to %s\n", outfilename_.c_str());
//  if (outfile_ != 0)
//    mprintf("\tOutput to file %s\n", outfile_->DataFilename().full());
//  if (!setname_.empty())
//    mprintf("\tData set name: %s\n", setname_.c_str());
//  masterDSL_ = DSL;

  return Action::OK;
}

// Action_CheckChirality::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_CheckChirality::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupCharMask(Mask1_)) return Action::ERR;
  if (Mask1_.None()) {
    mprinterr("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
    return Action::ERR;
  }
  // Reset any existing ResidueInfos to inactive
  for (Rarray::iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri)
    ri->isActive_ = 0;
  int Nactive = 0;
  // Loop over all non-solvent residues
  int resnum = 0;
  std::vector<std::string> NotFound;
  for (Topology::res_iterator res = currentParm->ResStart();
                              res != currentParm->ResEnd(); ++res, ++resnum)
  {
    int firstAtom = res->FirstAtom();
    int molnum = (*currentParm)[firstAtom].MolNum();
    // Skip solvent
    if (!currentParm->Mol(molnum).IsSolvent())
    {
      if (Mask1_.AtomsInCharMask( firstAtom, res->LastAtom() ))
      { 
        int n_atom = currentParm->FindAtomInResidue(resnum, "N");
        int ca_atom = currentParm->FindAtomInResidue(resnum, "CA");
        int c_atom = currentParm->FindAtomInResidue(resnum, "C");
        int cb_atom = currentParm->FindAtomInResidue(resnum, "CB");
        if (n_atom == -1 || ca_atom == -1 || c_atom == -1 || cb_atom == -1)
        {
          NotFound.push_back( currentParm->TruncResNameNum(resnum) );
          continue;
        }
        Nactive++;
        // See if a data set is already present
        Rarray::iterator rinfo = resInfo_.end();
        for (Rarray::iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri)
          if (resnum == ri->num_) {
            rinfo = ri;
            break;
          }
        if (rinfo != resInfo_.end()) {
          // Update coord indices in ResidueInfo
          rinfo->isActive_ = 1;
          rinfo->n_ = n_atom * 3;
          rinfo->ca_ = ca_atom * 3;
          rinfo->c_ = c_atom * 3;
          rinfo->cb_ = cb_atom * 3;
        } else {
          // New ResidueInfo
          ResidueInfo RI;
          RI.num_ = resnum;
          RI.isActive_ = 1;
          RI.n_ = n_atom * 3;
          RI.ca_ = ca_atom * 3;
          RI.c_ = c_atom * 3;
          RI.cb_ = cb_atom * 3;
          RI.N_L_ = 0;
          RI.N_D_ = 0;
          resInfo_.push_back( RI );
        }
      } // END atoms in mask
    } // END not solvent
  } // END loop over residues

  if (Nactive == 0) {
    mprintf("Warning: No valid residues selected from '%s'\n", currentParm->c_str());
    return Action::ERR;
  }
  mprintf("\tChecking chirality for %i residues\n", Nactive);
  if (!NotFound.empty()) {
    mprintf("\tSome atoms not found for %zu residues (this is expected for e.g. GLY)\n\t",
            NotFound.size());
    for (std::vector<std::string>::const_iterator rn = NotFound.begin();
                                                  rn != NotFound.end(); ++rn)
      mprintf(" %s", rn->c_str());
    mprintf("\n");
  }
  return Action::OK;
}

// Action_CheckChirality::DoAction()
Action::RetType Action_CheckChirality::DoAction(int frameNum, Frame* currentFrame,
                                            Frame** frameAddress)
{
  for (Rarray::iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri)
  {
    double torsion = Torsion( currentFrame->CRD(ri->n_),
                              currentFrame->CRD(ri->ca_),
                              currentFrame->CRD(ri->c_),
                              currentFrame->CRD(ri->cb_) );
    if (torsion < 0.0)
      ri->N_L_++;
    else
      ri->N_D_++;
  }

  return Action::OK;
}

void Action_CheckChirality::Print() {
  if (outfilename_.empty())
    mprintf("CHECKCHIRALITY: '%s'\n", Mask1_.MaskString());
  CpptrajFile out;
  if (out.OpenEnsembleWrite( outfilename_, ensembleNum_ )) return;
  out.Printf("%-8s %8s %8s\n", "#Res", "#L", "#D");
  for (Rarray::const_iterator ri = resInfo_.begin(); ri != resInfo_.end(); ++ri)
    out.Printf("%8i %8i %8i\n", ri->num_+1, ri->N_L_, ri->N_D_);
}
