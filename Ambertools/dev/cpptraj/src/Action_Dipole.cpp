#include <cmath> // sqrt
#include <algorithm> // std::max_element
#include "Action_Dipole.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL

// CONSTRUCTOR
Action_Dipole::Action_Dipole() :
  grid_(0),
  max_(0),
  CurrentParm_(0),
  ensembleNum_(-1)
{}

void Action_Dipole::Help() {
  mprintf("\t<filename>\n%s\n", GridAction::HelpText);
  mprintf("\t<mask1> {origin | box} [max <max_percent>]\n");
}

// Action_Dipole::Init()
Action::RetType Action_Dipole::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  // Get output filename
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: Dipole: no filename specified.\n");
    return Action::ERR;
  }
  // 'negative' means something different here than for other grid actions,
  // so get it here. Done this way to be consistent with PTRAJ behavior.
  if (actionArgs.hasKey("negative"))
    max_ = 1;
  else
    max_ = actionArgs.getKeyDouble("max", 0);
  // Get grid options
  grid_ = GridInit( "Dipole", actionArgs, *DSL );
  if (grid_ == 0) return Action::ERR;
  // Setup dipole x, y, and z grids
  dipole_.resize( grid_->Size(), Vec3(0.0,0.0,0.0) );

  // Get mask
  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: Dipole: No mask specified.\n");
    return Action::ERR;
  }
  mask_.SetMaskString(maskexpr);

  // Info
  mprintf("    DIPOLE:\n");
  GridInfo( *grid_ );
  mprintf("\tGrid will be printed to file %s\n",filename_.c_str());
  mprintf("\tMask expression: [%s]\n",mask_.MaskString());
  if (max_ > 0)
    mprintf("\tOnly keeping density >= to %.0lf%% of the maximum density\n", max_);

  return Action::OK;
}

// Action_Dipole::setup()
Action::RetType Action_Dipole::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->Nsolvent() < 1) {
    mprinterr("Error: Dipole: no solvent present in %s.\n", currentParm->c_str());
    return Action::ERR;
  }
  // Traverse over solvent molecules to find out the 
  // "largest" solvent molecule; allocate space for this
  // many coordinates.
  int NsolventAtoms = 0;
  for (Topology::mol_iterator Mol = currentParm->MolStart();
                              Mol != currentParm->MolEnd(); ++Mol)
  {
    if ( (*Mol).IsSolvent() ) {
      if ( (*Mol).NumAtoms() > NsolventAtoms )
        NsolventAtoms = (*Mol).NumAtoms();
    }
  }
  //sol_.resize( NsolventAtoms );
  mprintf("\tLargest solvent mol is %i atoms.\n", NsolventAtoms);

  // Setup grid, checks box info.
  if (GridSetup( *currentParm )) return Action::ERR;

  // Setup mask
  if (currentParm->SetupCharMask( mask_ ))
    return Action::ERR;
  mprintf("\t[%s] %i atoms selected.\n", mask_.MaskString(), mask_.Nselected());
  if (mask_.None()) {
    mprinterr("Error: Dipole: No atoms selected for parm %s\n", currentParm->c_str());
    return Action::ERR;
  }
  CurrentParm_ = currentParm;
  return Action::OK;
}

// Action_Dipole::action()
Action::RetType Action_Dipole::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Vec3 cXYZ, dipolar_vector, COM;

  // Set up center to origin or box center
  if (GridMode() == GridAction::BOX) 
    cXYZ = currentFrame->BoxCrd().Center();
  else if (GridMode() == GridAction::MASKCENTER)
    cXYZ = currentFrame->VGeometricCenter( CenterMask() );
  else // GridAction::ORIGIN/SPECIFIEDCENTER
    cXYZ.Zero();

  // Traverse over solvent molecules.
  //int i_solvent = 0; // DEBUG
  for (Topology::mol_iterator solvmol = CurrentParm_->MolStart();
                              solvmol != CurrentParm_->MolEnd(); ++solvmol)
  {
    if (!(*solvmol).IsSolvent()) continue;
    //++i_solvent; // DEBUG
    dipolar_vector.Zero();
    COM.Zero();
    double total_mass = 0;
    // Loop over solvent atoms
    for (int satom = (*solvmol).BeginAtom(); satom < (*solvmol).EndAtom(); ++satom)
    {
      if ( mask_.AtomInCharMask(satom) ) {
        // Get coordinates and shift to origin and then to appropriate spacing
        // NOTE: Do not shift into grid coords until the very end.
        const double* sol = currentFrame->XYZ( satom );
        // Calculate dipole vector. The center of mass of the solvent is used 
        // as the "origin" for the vector.
        // NOTE: the total charge on the solvent should be neutral for this 
        //       to have any meaning.
        double mass = (*CurrentParm_)[satom].Mass();
        total_mass += mass;
        COM[0] += (mass * sol[0]);
        COM[1] += (mass * sol[1]);
        COM[2] += (mass * sol[2]);

        double charge = (*CurrentParm_)[satom].Charge();
        dipolar_vector[0] += (charge * sol[0]);
        dipolar_vector[1] += (charge * sol[1]);
        dipolar_vector[2] += (charge * sol[2]);
      }
    }
    // If no atoms selected for this solvent molecule, skip.
    if (total_mass < Constants::SMALL) continue;

    // Grid COM
    COM /= total_mass;
    COM -= cXYZ;
    int ix, jy, kz;
    //mprintf("CDBG: Solvent %i XYZ %8.3f %8.3f %8.3f\n",solvmol-CurrentParm_->MolStart(),COM[0],COM[1],COM[2]);
    if (grid_->CalcBins( COM[0], COM[1], COM[2], ix, jy, kz )) {
      // Point COM is inside the grid. Increment grid and grid the dipole.
      long int bin = grid_->Increment( ix, jy, kz, Increment() );
      dipole_[bin] += dipolar_vector;
      //mprintf("CDBG: Indices %i %i %i\n", ix, jy, kz);
      //mprintf("CDBG: Bin = %lu\n", bin); 
    }
  } // END loop over solvent molecules

  return Action::OK;
}

// Action_Dipole::print()
/** Print dipole data in format for Chris Bayly's discern delegate that 
  * comes with Midas/Plus.
  */
void Action_Dipole::Print() {
  double max_density;
  CpptrajFile outfile;

  if (outfile.OpenEnsembleWrite(filename_, ensembleNum_)) {
    mprinterr("Error: Dipole: Cannot open output file.\n");
    return;
  }

  // Write header
  outfile.Printf("field 8\nsize 1\nnside 3\nnlayer 1\ndirectional\nvector\ndata\n");

  // Determine max density
  float maxF = *std::max_element(grid_->begin(), grid_->end());
  mprintf("\tDipole: maximum density is %f\n", maxF);

  if ( max_ > 0) {
    max_density = max_ * (double)maxF / 100.0;
    mprintf("\tWriting density if >= to %lf\n", max_density);
  } else
    max_density = 1.0;

  // Write data
  for (size_t k = 0; k < grid_->NZ(); ++k) {
    for (size_t j = 0; j < grid_->NY(); ++j) {
      for (size_t i = 0; i < grid_->NX(); ++i) {
        double density = grid_->GetElement(i, j, k);
        //mprintf("CDBG: %5i %5i %5i %lf\n",i,j,k,density);
        if ( density >= max_density ) {
          // Print Bin Coords
          Vec3 binCorner = grid_->BinCorner(i, j, k);
          outfile.Printf("%8.3f %8.3f %8.3f", binCorner[0], binCorner[1], binCorner[2]);
          // Normalize dipoles by density
          size_t idx = grid_->CalcIndex(i,j,k);
          dipole_[idx] /= density;
          // Write dipole components and length
          outfile.Printf(" %8.3f %8.3f %8.3f", 
                         dipole_[idx][0], dipole_[idx][1], dipole_[idx][2]);
          outfile.Printf(" %8.3f %8.3f\n", sqrt(dipole_[idx].Magnitude2()), density);
        }
      }
    }
  }
  outfile.CloseFile();
} 
