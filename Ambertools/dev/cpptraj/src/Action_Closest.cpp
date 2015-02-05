// Find closest waters to atoms in mask.
#include <cmath>
#include <algorithm> // sort
#include <cfloat> // DBL_MAX
#ifdef _OPENMP
#  include "omp.h"
#endif
#include "Action_Closest.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR
Action_Closest::Action_Closest() :
  outFile_(0),
  framedata_(0),
  moldata_(0),
  distdata_(0),
  atomdata_(0),
  Nclosest_(0),
  closestWaters_(0),
  firstAtom_(false),
  newParm_(0),
  NsolventMolecules_(0),
  debug_(0)
{} 

void Action_Closest::Help() {
  mprintf("\t<# to keep> <mask> [noimage] [first | oxygen]\n"
          "\t[closestout <filename> [name <setname>]] [outprefix <parmprefix>]\n"
          "  Keep only the closest <# to keep> solvent molecules to atoms in <mask>.\n"
          "  Molecules can be marked as solvent with the 'solvent' command.\n");
}

// DESTRUCTOR
Action_Closest::~Action_Closest() {
  if (newParm_!=0) delete newParm_;
}

// Action_Closest::Init()
Action::RetType Action_Closest::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  closestWaters_ = actionArgs.getNextInteger(-1);
  if (closestWaters_ < 0) {
    mprinterr("Error: Invalid # solvent molecules to keep (%i).\n",
              closestWaters_);
    return Action::ERR;
  }
  if ( actionArgs.hasKey("oxygen") || actionArgs.hasKey("first") )
    firstAtom_=true;
  InitImaging( !(actionArgs.hasKey("noimage")) );
  prefix_ = actionArgs.GetStringKey("outprefix");
  // Setup output file and sets if requested.
  // Will keep track of Frame, Mol#, Distance, and first solvent atom
  std::string filename = actionArgs.GetStringKey("closestout");
  if (!filename.empty()) {
    std::string dsetName = actionArgs.GetStringKey("name");
    if (dsetName.empty())
      dsetName = DSL->GenerateDefaultName("CLOSEST");
    // Set up datasets
    framedata_ = DSL->AddSetAspect(DataSet::INTEGER, dsetName, "Frame");
    moldata_   = DSL->AddSetAspect(DataSet::INTEGER, dsetName, "Mol");
    distdata_  = DSL->AddSetAspect(DataSet::DOUBLE,  dsetName, "Dist");
    atomdata_  = DSL->AddSetAspect(DataSet::INTEGER, dsetName, "FirstAtm");
    if (framedata_==0 || moldata_==0 || distdata_==0 || atomdata_==0) {
      mprinterr("Error: Could not setup data sets for output file %s\n",
                filename.c_str());
      return Action::ERR;
    }
    // Add sets to datafile in list.
    outFile_ = DFL->AddDataFile( filename );
    if (outFile_ == 0) {
      mprinterr("Error: Could not set up output file %s\n", filename.c_str());
      return Action::ERR;
    }
    outFile_->AddSet(framedata_);
    outFile_->AddSet(moldata_);
    outFile_->AddSet(distdata_);
    outFile_->AddSet(atomdata_);
    outFile_->ProcessArgs("noxcol");
  }

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty()) {
    mprinterr("Error: No mask specified.\n");
    return Action::ERR;
  }
  distanceMask_.SetMaskString(mask1);

  mprintf("    CLOSEST: Finding closest %i solvent molecules to atoms in mask %s\n",
          closestWaters_, distanceMask_.MaskString());
  if (!UseImage()) 
    mprintf("\tImaging will be turned off.\n");
  if (firstAtom_)
    mprintf("\tOnly first atom of solvent molecule used for distance calc.\n");
  if (outFile_!=0)
    mprintf("\tClosest molecules will be saved to %s\n",outFile_->DataFilename().base());
  if (!prefix_.empty())
    mprintf("\tStripped topology file will be written with prefix %s\n",
            prefix_.c_str());
  return Action::OK;
}

// Action_Closest::Setup()
/** Like the strip action, closest will modify the current parm keeping info
  * for atoms in mask plus the closestWaters solvent molecules. Set up the
  * vector of MolDist objects, one for every solvent molecule in the original
  * parm file. Atom masks for each solvent molecule will be set up.
  */
Action::RetType Action_Closest::Setup(Topology* currentParm, Topology** parmAddress) {
  // If there are no solvent molecules this action is not valid.
  if (currentParm->Nsolvent()==0) {
    mprintf("Warning: Parm %s does not contain solvent.\n",currentParm->c_str());
    return Action::ERR;
  }
  // If # solvent to keep >= solvent in this parm the action is not valid.
  if (closestWaters_ >= currentParm->Nsolvent()) {
    mprinterr("Error: # solvent to keep (%i) >= # solvent molecules in '%s' (%i)\n",
              closestWaters_, currentParm->c_str(), currentParm->Nsolvent());
    return Action::ERR;
  }
  SetupImaging( currentParm->BoxType() );
  if (ImagingEnabled())
    mprintf("\tDistances will be imaged.\n");
  else
    mprintf("\tImaging off.\n"); 
  // LOOP OVER MOLECULES
  // 1: Check that all solvent molecules contain same # atoms. Solvent 
  //    molecules must be identical for the command to work properly; 
  //    the prmtop strip occurs only once so the solvent params become fixed.
  // 2: Set up a mask for all solvent molecules.
  SolventMols_.clear();
  // NOTE: May not be necessary to init 'solvent'
  MolDist solvent;
  solvent.D = 0.0;
  solvent.mol = 0;
  SolventMols_.resize(currentParm->Nsolvent(), solvent);
  std::vector<MolDist>::iterator mdist = SolventMols_.begin();
  // 3: Set up the soluteMask for all non-solvent molecules.
  stripMask_.ResetMask();
  int molnum = 1;
  int newnatom = 0;
  int nclosest = 0;
  int NsolventAtoms = -1;
  keptWaterAtomNum_.resize(closestWaters_);
  for (Topology::mol_iterator Mol = currentParm->MolStart();
                              Mol != currentParm->MolEnd(); ++Mol)
  {
    if ( !Mol->IsSolvent() ) { // Not solvent, add to solute mask.
      stripMask_.AddAtomRange( Mol->BeginAtom(), Mol->EndAtom() );
      newnatom += Mol->NumAtoms();
    } else {                   // Solvent, check for same # of atoms.
      if (NsolventAtoms == -1)
        NsolventAtoms = Mol->NumAtoms();
      else if ( NsolventAtoms != Mol->NumAtoms() ) {
        mprinterr("Error: Solvent molecules in '%s' are not of uniform size.\n"
                  "Error:   First solvent mol = %i atoms, solvent mol %i = %i atoms.\n",
                  currentParm->c_str(), NsolventAtoms, molnum, (*Mol).NumAtoms());
        return Action::ERR;
      }
      // mol here is the output molecule number which is why it starts from 1.
      mdist->mol = molnum;
      // Solvent molecule mask
      mdist->mask.AddAtomRange( Mol->BeginAtom(), Mol->EndAtom() );
      // Atoms in the solvent molecule to actually calculate distances to.
      if (firstAtom_) {
        mdist->solventAtoms.assign(1, Mol->BeginAtom() );
      } else {
        mdist->solventAtoms.clear();
        mdist->solventAtoms.reserve( Mol->NumAtoms() );
        for (int svatom = Mol->BeginAtom(); svatom < Mol->EndAtom(); svatom++)
          mdist->solventAtoms.push_back( svatom );
      }
      // For solvent molecules that will be kept, record what the atom number
      // will be in the new stripped parm.
      if (nclosest < closestWaters_) {
        keptWaterAtomNum_[nclosest] = newnatom;
        stripMask_.AddAtomRange( Mol->BeginAtom(), Mol->EndAtom() );
        newnatom += Mol->NumAtoms();
        ++nclosest;
      }
      //SolventMols[solventMol].mask.PrintMaskAtoms("solvent");
      ++mdist;
    }
    ++molnum;
  }

  // Setup distance atom mask
  // NOTE: Should ensure that no solvent atoms are selected!
  if ( currentParm->SetupIntegerMask(distanceMask_) ) return Action::ERR;
  if (distanceMask_.None()) {
    mprintf("Warning: Distance mask '%s' contains no atoms.\n",
            distanceMask_.MaskString());
    return Action::ERR;
  }
  distanceMask_.MaskInfo();

  // Check the total number of solvent atoms to be kept.
  NsolventAtoms *= closestWaters_;
  mprintf("\tKeeping %i solvent atoms.\n",NsolventAtoms);
  if (NsolventAtoms < 1) {
    mprinterr("Error: # of solvent atoms to be kept is < 1.\n");
    return Action::ERR;
  }
  // Store original # of molecules.
  // NOTE: This is stored so that it can be used in the OpenMP section
  //       of action. I dont think iterators are thread-safe.
  NsolventMolecules_ = currentParm->Nsolvent();
  // Create stripped Parm
  if (newParm_!=0) delete newParm_;
  newParm_ = currentParm->modifyStateByMask(stripMask_);
  if (newParm_==0) {
    mprinterr("Error: Could not create new parmtop.\n");
    return Action::ERR;
  }
  newParm_->Brief("Closest parm:");
  // Allocate space for new frame
  newFrame_.SetupFrameV( newParm_->Atoms(), newParm_->HasVelInfo(), newParm_->NrepDim() );
  // If prefix given then output stripped parm
  if (!prefix_.empty()) {
    ParmFile pfile;
    if ( pfile.WritePrefixTopology(*newParm_, prefix_, ParmFile::AMBERPARM, debug_ ) )
      mprinterr("Error: Could not write out 'closest' parm file.\n");
  }
  // Set parm
  *parmAddress = newParm_;

  return Action::OK;  
}

// Action_Closest::DoAction()
/** Find the minimum distance between atoms in distanceMask and each 
  * solvent Mask.
  */
Action::RetType Action_Closest::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress)
{
  int solventMol; 
  double Dist, maxD;
  Matrix_3x3 ucell, recip;
  AtomMask::const_iterator solute_atom;
  Iarray::const_iterator solvent_atom;

  if (ImagingEnabled()) {
    currentFrame->BoxCrd().ToRecip(ucell, recip);
    // Calculate max possible imaged distance
    maxD = currentFrame->BoxCrd().BoxX() + currentFrame->BoxCrd().BoxY() + 
           currentFrame->BoxCrd().BoxZ();
    maxD *= maxD;
  } else {
    // If not imaging, set max distance to an arbitrarily large number
    maxD = DBL_MAX;
  }

  // Loop over all solvent molecules in original frame
  // DEBUG
  //mprintf("Closest: Begin parallel loop for %i\n",frameNum);
#ifdef _OPENMP
#pragma omp parallel private(solventMol,solute_atom,Dist,solvent_atom)
{
  //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
#pragma omp for
#endif
  for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
    //mprintf("[%i] Calculating distance for molecule %i\n",omp_get_thread_num(),solventMol);
    // Set the initial minimum distance for this solvent mol to be the
    // max possible distance.
    SolventMols_[solventMol].D = maxD;
    // Calculate distance between each atom in distanceMask and atoms in solvent Mask
    for (solvent_atom = SolventMols_[solventMol].solventAtoms.begin();
         solvent_atom != SolventMols_[solventMol].solventAtoms.end(); ++solvent_atom)
    {
      for (solute_atom = distanceMask_.begin(); 
           solute_atom != distanceMask_.end(); ++solute_atom)
      {
        Dist = DIST2(currentFrame->XYZ(*solute_atom),
                     currentFrame->XYZ(*solvent_atom), ImageType(), 
                     currentFrame->BoxCrd(), ucell, recip);
        if (Dist < SolventMols_[solventMol].D) 
          SolventMols_[solventMol].D = Dist;
        //mprintf("DBG: SolvMol %i, soluteAtom %i, solventAtom %i, D= %f, minD= %f\n",
        //        solventMol, *solute_atom, *solvent_atom, Dist, sqrt(SolventMols_[solventMol].D));
      }
    }
    // DEBUG - Print distances
    //mprintf("DEBUG:\tMol %8i minD= %lf\n",solventMol, SolventMols[solventMol].D);
  } // END for loop over solventMol
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  // DEBUG
  //mprintf("Closest: End parallel loop for %i, got %i Distances.\n",frameNum,(int)SolventMols.size());
  // DEBUG

  // Sort distances
  std::sort( SolventMols_.begin(), SolventMols_.end(), moldist_cmp() );
  // Add first closestWaters solvent atoms to stripMask
  std::vector<MolDist>::iterator solventend = SolventMols_.begin() + closestWaters_;
  Iarray::const_iterator katom = keptWaterAtomNum_.begin();
  for ( std::vector<MolDist>::const_iterator solvent = SolventMols_.begin();
                                             solvent != solventend;
                                           ++solvent, ++katom ) 
  {
    //mprintf("DEBUG:\tmol %i ",(*solvent).mol);
    //solvent->mask.PrintMaskAtoms("Mask");
    stripMask_.AddMaskAtPosition( solvent->mask, *katom );
    // Record which water molecules are closest if requested
    if (outFile_!=0) {
      int fnum = frameNum + 1;
      framedata_->Add(Nclosest_, &fnum);
      moldata_->Add(Nclosest_, &(solvent->mol));
      Dist = sqrt( solvent->D );
      distdata_->Add(Nclosest_, &Dist);
      solvent_atom = solvent->mask.begin();
      int solvent_first_atom = *solvent_atom + 1; 
      atomdata_->Add(Nclosest_, &solvent_first_atom);
      ++Nclosest_;
    }
    // DEBUG - print first closestWaters distances
    //mprintf("DEBUG: Mol %i   D2= %lf   Atom0= %i\n",(*it).mol, (*it).D, (*it).mask->Selected[0]);
  }

  // Modify and set frame
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("action_stripMask");
  newFrame_.SetFrame(*currentFrame, stripMask_);
  *frameAddress = &newFrame_;

  return Action::OK;
}
