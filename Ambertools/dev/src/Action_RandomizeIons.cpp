#include <cmath> // sqrt
#include <cstdlib> // random, srandom
#include "Action_RandomizeIons.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_RandomizeIons::Action_RandomizeIons() :
  overlap_(0),
  min_(0),
  seed_(1),
  debug_(0),
  n_solvent_(0)
{}

void Action_RandomizeIons::Help() {
  mprintf("\t<mask> [around <mask> by <distance>] [overlap <value>]\n"
          "\t[noimage] [seed <value>]\n"
          "  Swap positions of ions in <mask> with randomly selected solvent molecules.\n");
}

// Action_RandomizeIons::Init()
Action::RetType Action_RandomizeIons::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get first mask
  std::string ionmask = actionArgs.GetMaskNext();
  if (ionmask.empty()) {
    mprinterr("Error: randomizeions: No mask for ions specified.\n");
    return Action::ERR;
  }
  ions_.SetMaskString( ionmask );

  // Get Keywords
  InitImaging( !actionArgs.hasKey("noimage") );
  seed_ = actionArgs.getKeyInt("seed", -1);
  overlap_ = actionArgs.getKeyDouble("overlap", 3.5);
  min_ = actionArgs.getKeyDouble("by", 3.5);
  // Pre-square overlap and min
  overlap_ *= overlap_;
  min_ *= min_;
  // If no around mask specified, leave blank
  aroundmask_ = actionArgs.GetStringKey("around");
  if (!aroundmask_.empty())
    around_.SetMaskString( aroundmask_ );
  
  // INFO
  mprintf("    RANDOMIZEIONS: swapping the postions of the ions in mask [%s]\n", 
          ions_.MaskString());
  mprintf("                   with the solvent. No ions can get closer than %.2f angstroms\n",
          sqrt( overlap_ ));
  mprintf("                   to another ion.\n");
  if (!aroundmask_.empty()) {
    mprintf("                   No ion can get closer than %.2f angstroms to mask [%s]\n",
            sqrt( min_ ), around_.MaskString());
  }
  if (!UseImage())
    mprintf("                   Imaging of the coordinates will not be performed.\n");
  if (seed_ > 0) {
    mprintf("                   Random number generator seed is %i\n", seed_);
    srandom((unsigned) seed_);
  }

  return Action::OK;
}

// Action_RandomizeIons::setup()
Action::RetType Action_RandomizeIons::Setup(Topology* currentParm, Topology** parmAddress) {
  n_solvent_ = currentParm->Nsolvent();
  if (n_solvent_==0) {
    mprinterr("Warning: randomizeions: This command only works if solvent information\n");
    mprinterr("Warning:                has been specified.");
    //mprinterr(" See the \"solvent\" command.");
    mprinterr("\n");
    return Action::ERR;
  }

  // Set up ion mask
  if (currentParm->SetupIntegerMask( ions_ )) return Action::ERR;
  if ( ions_.None() ) {
    mprinterr("Warning: randomizeions: Mask [%s] has no atoms.\n", ions_.MaskString());
    return Action::ERR;
  }
  mprintf("\tIon mask is [%s] (%i atoms)\n", ions_.MaskString(), ions_.Nselected());

  // Set up the around mask if necessary
  if (!aroundmask_.empty()) {
    if (currentParm->SetupIntegerMask( around_ )) return Action::ERR;
    if ( around_.None() ) {
      mprintf("Warning: randomizeions: Around mask [%s] has no atoms.\n",
              around_.MaskString());
    } else {
      mprintf("\tAround mask is [%s] (%i atoms)\n", around_.MaskString(),
              around_.Nselected());
    }
  }

  // Check that each ion is only a single atom residue.
  // NOTE: Should this be a molecule check instead? If so can then get rid of ResSize 
  for (AtomMask::const_iterator ion = ions_.begin(); ion != ions_.end(); ++ion)
  {
    int res = (*currentParm)[*ion].ResNum();
    if (debug_ > 0)
      mprintf("\tAtom %i is in residue %i which is %i atoms\n",
              *ion+1, res+1, currentParm->Res( res ).NumAtoms() );
    if ( currentParm->Res( res ).NumAtoms() > 1 ) {
      mprintf("Warning: randomizeions: Ion atom %i belongs to residue %i which\n",
              *ion + 1, res + 1);
      mprintf("Warning:                contains more than 1 atom (%i)!\n", 
              currentParm->Res( res ).NumAtoms());
    }
  }

  // Check the solvent information to make sure that each solvent listed has the
  // same number of atoms in each molecule; otherwise a uniform trajectory is not
  // possible and therefore this command will be ignored.
  // Also save the start and end atom of each solvent molecule.
  int NsolventAtoms = -1;
  solventStart_.clear();
  solventEnd_.clear();
  solventStart_.reserve( n_solvent_ );
  solventEnd_.reserve( n_solvent_ );
  for (Topology::mol_iterator Mol = currentParm->MolStart();
                              Mol != currentParm->MolEnd(); ++Mol)
  {
    if ( (*Mol).IsSolvent() ) {
      if (NsolventAtoms == -1)
        NsolventAtoms = (*Mol).NumAtoms();
      else if ( NsolventAtoms != (*Mol).NumAtoms() ) {

        mprinterr("Warning: randomizeions: Solvent molecules in %s are not of uniform size.\n",
                  currentParm->c_str());
        mprinterr("       First solvent mol = %i atoms, this solvent mol = %i atoms.\n",
                  NsolventAtoms, (*Mol).NumAtoms());
        return Action::ERR;
      }
      solventStart_.push_back( (*Mol).BeginAtom() );
      solventEnd_.push_back( (*Mol).EndAtom() );
    }
  }
  SetupImaging( currentParm->BoxType() );
  // Allocate solvent molecule mask
  solvent_.resize( n_solvent_ );

  return Action::OK;
}

// Action_RandomizeIons::action()
Action::RetType Action_RandomizeIons::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 ucell, recip;
  Vec3 trans;

  if (ImageType()==NONORTHO)
    currentFrame->BoxCrd().ToRecip(ucell, recip);
  // loop over all solvent molecules and mark those that are too close to the solute
  std::vector<bool>::iterator smask = solvent_.begin();
  //int smolnum = 1; // DEBUG
  for (std::vector<int>::iterator beginatom = solventStart_.begin();
                                  beginatom != solventStart_.end(); ++beginatom)
  {
    *smask = true;
    // is solvent molecule to near any atom in the around mask?
    if (!aroundmask_.empty()) {
      for (AtomMask::const_iterator atom = around_.begin(); atom != around_.end(); ++atom)
      {
        double dist = DIST2( currentFrame->XYZ(*beginatom), currentFrame->XYZ(*atom), 
                             ImageType(), currentFrame->BoxCrd(), ucell, recip);
        //mprintf("CDBG: @%i to @%i = %lf\n", *beginatom+1,
        //        *atom+1, sqrt(dist));
        if (dist < min_) {
          *smask = false;
          //mprintf("RANDOMIZEIONS: water %i only %.2f ang from around @%i\n",
          //        smolnum, sqrt(dist), *atom+1);
          break;
        }
      }
    }
    ++smask;
    //++smolnum; // DEBUG
  }

  // DEBUG - print solvent molecule mask
  //mprintf("RANDOMIZEIONS: The following waters are ACTIVE so far:\n");
  //smolnum = 1; // DEBUG
  int smoltot = 0;
  for (smask = solvent_.begin(); smask != solvent_.end(); ++smask) {
    if (*smask) {
      //mprintf(" %5i ", smolnum);
      ++smoltot;
      //if (smoltot%10 == 0) mprintf("\n");
    }
    //++smolnum; // DEBUG
  }
  if (debug_ > 2)
    mprintf("RANDOMIZEIONS: A total of %i waters (out of %zu) are active\n",
            smoltot, solvent_.size());

  // Outer loop over all ions
  for (AtomMask::const_iterator ion = ions_.begin(); ion != ions_.end(); ++ion)
  {
    //mprintf("RANDOMIZEIONS: Processing ion atom %i\n", *ion+1);
    // is a potential solvent molecule close to any of the ions (except this one)?
    smask = solvent_.begin();
    //smolnum = 1; // DEBUG
    for (std::vector<int>::iterator beginatom = solventStart_.begin();
                                    beginatom != solventStart_.end(); ++beginatom)
    {
      if (*smask) {
        // if this solvent is active, check distance to all other ions
        for (AtomMask::const_iterator ion2 = ions_.begin(); ion2 != ions_.end(); ++ion2)
        {
          if (*ion != *ion2) {
            double dist = DIST2( currentFrame->XYZ(*beginatom), currentFrame->XYZ(*ion2), 
                                 ImageType(), currentFrame->BoxCrd(), ucell, recip);
            if (dist < overlap_) {
              *smask = false;
              //mprintf("RANDOMIZEIONS: water %i only %.2f ang from ion @%i\n",
              //        smolnum, sqrt(dist), *ion2+1);
              break;
            }
          }
        } // END inner loop over ions
      }
      ++smask;
      //++smolnum; // DEBUG
    } // END loop over solvent molecules

    // solvent should now be true for all solvent molecules eligible to
    // swap with ion.
    int loop = 1;
    int swapMol = 0;
    while (loop > 0 && loop < 10000) {
      // Run the random number generator so that the same number is not produced
      // when the seed was set manually.
      swapMol = random() % n_solvent_;
      if ( solvent_[swapMol] ) 
        loop = -1;
      else 
        ++loop;
    }

    // If a suitable solvent molecule was found, swap it.
    if (loop > 0) {
      mprintf("Warning: randomizeions: Tried to swap ion @%i with %i random waters\n",*ion+1,loop);
      mprintf("Warning:                and couldn't meet criteria; skipping.\n");
    } else {
      if (debug_ > 2)
        mprintf("RANDOMIZEIONS: Swapping solvent mol %i for ion @%i\n", swapMol+1, *ion+1);
      const double* ionXYZ = currentFrame->XYZ( *ion );
      int sbegin = solventStart_[ swapMol ];
      const double* watXYZ = currentFrame->XYZ( sbegin );
      Vec3 trans( ionXYZ[0]-watXYZ[0], ionXYZ[1]-watXYZ[1], ionXYZ[2]-watXYZ[2]);
      currentFrame->Translate( trans, sbegin, solventEnd_[ swapMol ] );
      trans.Neg();
      currentFrame->Translate( trans, *ion );
    }

  } // END outer loop over all ions

  return Action::OK;
}
