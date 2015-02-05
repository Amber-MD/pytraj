#include <cmath> // exp
#include <algorithm> // std::max_element
#include "Action_Grid.h"
#include "CpptrajStdio.h"
#include "PDBfile.h"

// CONSTRUCTOR
Action_Grid::Action_Grid() :
  normalize_(NONE),
  ensembleNum_(-1),
  // Default particle density (molecules/Ang^3) for water based on 1.0 g/mL
  density_(0.033456),
  max_(0.80),
  madura_(0),
  smooth_(0),
  nframes_(0),
  invert_(false),
  grid_(0)
{}

void Action_Grid::Help() {
  mprintf("\t<filename>\n%s\n", GridAction::HelpText);
  mprintf("\t<mask> [[smoothdensity <value>] [invert]] [madura <madura>]\n"
          "\t[pdb <pdbout> [max <fraction>]] [normframe | normdensity [density <density>]]\n"
          "  Bin atoms in <mask> into a 3D grid.\n"
          "    <fraction>: Percent of max to write.\n"
          "    <madura>  : Grid values lower than <madura> become flipped in sign, exposes low density.\n"
          "    <value>   : Used to smooth density.\n");
}

// Action_Grid::Init()
Action::RetType Action_Grid::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  nframes_ = 0;
  // Get output filename
  std::string filename = actionArgs.GetStringNext();
  if (filename.empty()) {
    mprinterr("Error: GRID: no filename specified.\n");
    return Action::ERR;
  }
  // Get grid options
  grid_ = GridInit( "GRID", actionArgs, *DSL );
  if (grid_ == 0) return Action::ERR;

  // Get extra options
  max_ = actionArgs.getKeyDouble("max", 0.80);
  madura_ = actionArgs.getKeyDouble("madura", 0);
  smooth_ = actionArgs.getKeyDouble("smoothdensity", 0);
  invert_ = actionArgs.hasKey("invert");
  pdbname_ = actionArgs.GetStringKey("pdb");
  density_ = actionArgs.getKeyDouble("density",0.033456);
  if (actionArgs.hasKey("normframe")) normalize_ = TO_FRAME;
  else if (actionArgs.hasKey("normdensity")) normalize_ = TO_DENSITY;
  else normalize_ = NONE;
  if (normalize_ != NONE && (smooth_ > 0.0 || madura_ > 0.0)) {
    mprinterr("Error: Normalize options are not compatible with smoothdensity/madura options.\n");
    return Action::ERR;
  }
  // Get mask
  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: GRID: No mask specified.\n");
    return Action::ERR;
  }
  mask_.SetMaskString(maskexpr);

  // Setup output file
  DataFile* outfile = DFL->AddDataFile(filename, actionArgs);
  if (outfile == 0) {
    mprinterr("Error: grid: Could not set up output file %s\n", filename.c_str());
    return Action::ERR;
  }
  outfile->AddSet((DataSet*)grid_);
  // grid_.PrintXplor( filename_, "This line is ignored", 
  //                      "rdparm generated grid density" );

  // Info
  mprintf("    GRID:\n");
  GridInfo( *grid_ );
  mprintf("\tGrid will be printed to file %s\n",filename.c_str());
  mprintf("\tMask expression: [%s]\n",mask_.MaskString());
  if (pdbname_.empty())
    mprintf("\tPseudo-PDB will be printed to STDOUT.\n");
  else
    mprintf("\tPseudo-PDB will be printed to %s\n", pdbname_.c_str());
  if (normalize_ == TO_FRAME)
    mprintf("\tGrid will be normalized by number of frames.\n");
  else if (normalize_ == TO_DENSITY)
    mprintf("\tGrid will be normalized to a density of %g molecules/Ang^3.\n", density_);
  // TODO: print extra options

  return Action::OK;
}

// Action_Grid::Setup()
Action::RetType Action_Grid::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup grid, checks box info.
  if (GridSetup( *currentParm )) return Action::ERR;

  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ ))
    return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: GRID: No atoms selected for parm %s\n", currentParm->c_str());
    return Action::ERR;
  }

  return Action::OK;
}

// Action_Grid::DoAction()
Action::RetType Action_Grid::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  GridFrame( *currentFrame, mask_, *grid_ );
  ++nframes_;
  return Action::OK;
}

// Action_Grid::print()
void Action_Grid::Print() {
  if (nframes_ < 1) return;
  // Perform normalization and find max.
  double gridMax = 0.0;
  if (normalize_ == NONE) {
    mprintf("    GRID: No normalization");
    if (smooth_ > 0.0) mprintf(", smoothing factor (%g)", smooth_);
    mprintf(".\n");
    for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval) {
      double gridval = (double)(*gval);
      // ----- SMOOTHING -----
      if (smooth_ > 0.0) {
        double yy = gridval - smooth_;
        double xx = yy*yy / (0.2 * smooth_ * smooth_);
        xx = exp( -xx );
        if (invert_) {
          if (gridval > smooth_) // NOTE: Comparison OK? Needs cast?
            gridval = -5.0;
          else
            gridval -= gridval * xx;
          /* COMMENTED OUT IN ORIGINAL PTRAJ CODE
          if (gridInfo->grid[index] < action->darg3) {
            gridInfo->grid[index] = 0.0;
          }
          */
          if (gridval >= 0)
            gridval = smooth_ - gridval;
        } else {
          if (gridval < smooth_)
            gridval = 0;
          else
            gridval -= gridval * xx;
          if (gridval < smooth_)
            gridval = 0;
        }
      }
      // do the madura negative option to expose low density
      if ( madura_ > 0.0 && gridval > 0.0 && gridval < madura_ )
        *gval = (float)-gridval;
      else
        *gval = (float) gridval;
      if ( gridval > gridMax )
        gridMax = gridval;
    }
  } else {
    // Normalize to frames / density.
    mprintf("    GRID: Normalization");
    double dens = 1.0;
    if (normalize_ == TO_DENSITY) {
      dens = grid_->VoxelVolume() * density_;
      mprintf(" to density %g mol/Ang^3, voxel volume= %g Ang^3, %g mols/voxel,",
              density_, grid_->VoxelVolume(), dens);
    } else
      mprintf(" to");
    mprintf(" number of frames %u", nframes_);
    double norm = 1.0 / ((double)nframes_ * dens);
    mprintf(", normalization factor= %g\n",norm);
    for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval) {
      double gridval = (double)(*gval) * norm;
      gridMax = std::max(gridval, gridMax);
      *gval = (float)gridval;
    }
  }
  mprintf("\tGrid max is %g\n", gridMax);
  // PDBfile output
  PrintPDB( gridMax );
}

// Action_Grid::PrintPDB()
void Action_Grid::PrintPDB(double gridMax)
{
  if (gridMax == 0.0) {
    mprinterr("Error: Grid max is 0. No density for PDB write.\n");
    return;
  }
  double norm = 1.0 / gridMax;
  // Calculate normalization if necessary
//  if (norm < 0.0) {
//    norm = (double)*std::max_element(grid_->begin(), grid_->end());
//  }
  // Write PDB
  PDBfile pdbout;
  if (pdbout.OpenEnsembleWrite(pdbname_, ensembleNum_)) {
    mprinterr("Error: Cannot open PDB for grid output.\n");
    return;
  }
  mprintf("\tWriting PDB of grid points > %.2f%% of grid max.\n", max_*100.0);
  int res = 1;
  for (size_t k = 0; k < grid_->NZ(); ++k) {
    for (size_t j = 0; j < grid_->NY(); ++j) {
      for (size_t i = 0; i < grid_->NX(); ++i) {
        double gridval = grid_->GetElement(i, j, k) * norm;
        if (gridval > max_) {
          Vec3 cxyz = grid_->BinCenter(i,j,k);
          pdbout.WriteATOM(res++, cxyz[0], cxyz[1], cxyz[2], "GRID", gridval);
        }
      }
    }
  }
  // Write grid boundaries
  for (size_t k = 0; k <= grid_->NZ(); k += grid_->NZ())
    for (size_t j = 0; j <= grid_->NY(); j += grid_->NY())
      for (size_t i = 0; i <= grid_->NX(); i += grid_->NX()) {
        Vec3 cxyz = grid_->BinCenter(i,j,k);
        pdbout.WriteHET(res, cxyz[0], cxyz[1], cxyz[2]);
      }
}
