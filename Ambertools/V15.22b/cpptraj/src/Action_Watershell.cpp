#include "Action_Watershell.h"
#include "CpptrajStdio.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_Watershell::Action_Watershell() :
  lowerCutoff_(0),
  upperCutoff_(0),
  CurrentParm_(0),
  lower_(0),
  upper_(0),
  numthreads_(1)
# ifdef _OPENMP
  ,activeResidues_thread_(0),
  NactiveResidues_(0)
# endif
{ }
#ifdef _OPENMP
Action_Watershell::~Action_Watershell() {
  if (activeResidues_thread_ != 0) {
    for (int i = 0; i < numthreads_; ++i)
      delete[] activeResidues_thread_[i];
    delete activeResidues_thread_;
  }
}
#endif

void Action_Watershell::Help() {
  mprintf("\t<solutemask> <filename> [lower <lower cut>] [upper <upper cut>]\n"
          "\t[noimage] [<solventmask>]\n"
          "  Calculate # of waters in 1st and 2nd solvation shells (defined by\n"
          "  <lower cut> (default 3.4 Ang.) and <upper cut> (default 5.0 Ang.)\n"
          "  distance cut-offs respectively.\n");
}

// Action_Watershell::init()
Action::RetType Action_Watershell::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  InitImaging( !actionArgs.hasKey("noimage") );

  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: WATERSHELL: Solute mask must be specified.\n");
    return Action::ERR;
  }
  soluteMask_.SetMaskString( maskexpr );

  std::string filename = actionArgs.GetStringNext();
  if (filename.empty()) {
    mprinterr("Error: WATERSHELL: Output filename must be specified.\n");
    return Action::ERR;
  }

  lowerCutoff_ = actionArgs.getKeyDouble("lower", 3.4);
  upperCutoff_ = actionArgs.getKeyDouble("upper", 5.0);

  // Check for solvent mask
  solventmaskexpr_ = actionArgs.GetMaskNext();

  // Set up datasets
  std::string dsname = actionArgs.GetStringNext();
  if (dsname.empty())
    dsname = DSL->GenerateDefaultName("WS");
  lower_ = DSL->AddSetAspect(DataSet::INTEGER, dsname, "lower");
  upper_ = DSL->AddSetAspect(DataSet::INTEGER, dsname, "upper");
  if (lower_ == 0 || upper_ == 0) return Action::ERR;
  DFL->AddSetToFile(filename, lower_);
  DFL->AddSetToFile(filename, upper_);
# ifdef _OPENMP
  // Determine number of parallel threads
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    numthreads_ = omp_get_num_threads();
}
# endif
  mprintf("    WATERSHELL: Output to %s\n",filename.c_str());
  if (!UseImage())
    mprintf("                Imaging is disabled.\n");
  mprintf("                The first shell will contain water < %.3lf angstroms from\n",
          lowerCutoff_);
  mprintf("                the solute; the second shell < %.3lf angstroms...\n",
          upperCutoff_);
  mprintf("                Solute atoms will be specified by [%s]\n",soluteMask_.MaskString());
  if (!solventmaskexpr_.empty()) {
    mprintf("                Solvent atoms will be specified by [%s]\n",
            solventmaskexpr_.c_str());
    solventMask_.SetMaskString( solventmaskexpr_ );
  }
  if (numthreads_ > 1)
    mprintf("\tParallelizing calculation with %i threads.\n", numthreads_);

  // Pre-square upper and lower cutoffs
  lowerCutoff_ *= lowerCutoff_;
  upperCutoff_ *= upperCutoff_;

  return Action::OK;
}

// Action_Watershell::setup()
/** Set up solute and solvent masks. If no solvent mask was specified use 
  * solvent information in the current topology.
  */
Action::RetType Action_Watershell::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up solute mask
  if (currentParm->SetupIntegerMask( soluteMask_ )) return Action::ERR;
  if ( soluteMask_.None() ) {
    mprinterr("Error: WATERSHELL: No atoms in solute mask [%s].\n",soluteMask_.MaskString());
    return Action::ERR;
  }
  // Set up solvent mask
  if (!solventmaskexpr_.empty()) {
    if (currentParm->SetupIntegerMask( solventMask_ )) return Action::ERR;
  } else {
    solventMask_.ResetMask();
    for (Topology::mol_iterator mol = currentParm->MolStart();
                                mol != currentParm->MolEnd(); ++mol)
    {
      if ( (*mol).IsSolvent() )
        solventMask_.AddAtomRange( (*mol).BeginAtom(), (*mol).EndAtom() );
    }
  }
  if ( solventMask_.None() ) {
    if (!solventmaskexpr_.empty())
      mprinterr("Error: WATERSHELL: No solvent atoms selected by mask [%s]\n",
                solventmaskexpr_.c_str());
    else
      mprinterr("Error: WATERSHELL: No solvent atoms in topology %s\n",currentParm->c_str());
    return Action::ERR;
  }
  SetupImaging( currentParm->BoxType() );
  // Create space for residues
# ifdef _OPENMP
  // Only re-allocate for larger # of residues
  if ( currentParm->Nres() > NactiveResidues_ ) {
    if (activeResidues_thread_ != 0) {
      // Deallocate each thread
      for (int i = 0; i < NactiveResidues_; ++i)
        delete[] activeResidues_thread_[i];
    } else {
      // Initial thread allocation needed
      activeResidues_thread_ = new int*[ numthreads_ ];
    }
    // Allocate each thread
    for (int i = 0; i < numthreads_; ++i) {
      activeResidues_thread_[i] = new int[ currentParm->Nres() ];
      std::fill( activeResidues_thread_[i], activeResidues_thread_[i] + currentParm->Nres(), 0 );
    }
  }
  NactiveResidues_ = currentParm->Nres();
# else
  activeResidues_.resize( currentParm->Nres(), 0 );
# endif
  // Store current Parm
  CurrentParm_ = currentParm;
  return Action::OK;    
}

// Action_Watershell::action()
Action::RetType Action_Watershell::DoAction(int frameNum, Frame* currentFrame, 
                                            Frame** frameAddress) 
{
  Matrix_3x3 ucell, recip;
  int nlower = 0;
  int nupper = 0;

  if (ImageType()==NONORTHO) currentFrame->BoxCrd().ToRecip(ucell,recip);

  int Vidx, Uidx, Vat, currentRes;
  int NU = soluteMask_.Nselected();
  int NV = solventMask_.Nselected();
  double dist;
# ifdef _OPENMP
  int mythread;
#pragma omp parallel private(dist,Vidx,Vat,currentRes,Uidx,mythread)
{
  mythread = omp_get_thread_num();
#pragma omp for
# endif
  // Assume solvent mask is the larger one.
  // Loop over solvent atoms
  for (Vidx = 0; Vidx < NV; Vidx++) {
    // Figure out which solvent residue this is
    Vat = solventMask_[Vidx];
    currentRes = (*CurrentParm_)[ Vat ].ResNum();
    // Loop over solute atoms
    for (Uidx = 0; Uidx < NU; Uidx++) {
      // If residue is not yet marked as 1st shell, calc distance
#     ifdef _OPENMP
      if ( activeResidues_thread_[mythread][currentRes] < 2 )
#     else
      if ( activeResidues_[currentRes] < 2 )
#     endif
      {
        dist = DIST2(currentFrame->XYZ(soluteMask_[Uidx]), currentFrame->XYZ(Vat),
                     ImageType(), currentFrame->BoxCrd(), ucell, recip );
        // Less than upper, 2nd shell
        if (dist < upperCutoff_) 
        {
#         ifdef _OPENMP
          activeResidues_thread_[mythread][currentRes] = 1;
#         else
          activeResidues_[currentRes] = 1;
#         endif
          // Less than lower, 1st shell
          if (dist < lowerCutoff_)
#           ifdef _OPENMP
            activeResidues_thread_[mythread][currentRes] = 2;
#           else
            activeResidues_[currentRes] = 2;
#           endif
        }
      }
    } // End loop over solute atoms
  } // End loop over solvent atoms
# ifdef _OPENMP
} // END parallel
  // Combine results from each thread.
  for (int res = 0; res < NactiveResidues_; res++) {
    int shell = 0;
    for (int thread = 0; thread < numthreads_; thread++) {
      if (activeResidues_thread_[thread][res] > shell)
        shell = activeResidues_thread_[thread][res];
      // Dont break here so we can reset counts. Could also do with a fill 
      activeResidues_thread_[thread][res] = 0;
    }
    if (shell > 0) {
      ++nupper;
      if (shell > 1) ++nlower;
    }
  }
# else
  // Now each residue is marked 0 (no shell), 1 (second shell), 2 (first shell)
  for (std::vector<int>::iterator shell = activeResidues_.begin();
                                  shell != activeResidues_.end(); ++shell)
  {
    if ( *shell > 0 ) {
      ++nupper;
      if ( *shell > 1 ) ++nlower;
    }
    // Reset for next pass
    *shell = 0;
  }
# endif
  lower_->Add(frameNum, &nlower);
  upper_->Add(frameNum, &nupper);

  return Action::OK;
}
