#include <cfloat>
#include <cmath> // ceil
#include "Action_Bounds.h"
#include "CpptrajStdio.h"
#include "DataSet_3D.h" // For allocating grid

// CONSTRUCTOR
Action_Bounds::Action_Bounds() : ensembleNum_(-1), offset_(1), grid_(0) {}

void Action_Bounds::Help() {
  mprintf("\t[<mask>] [out <filename>] [dx <dx> [dy <dy>] [dz <dz>] name <gridname>]\n"
          "\t[offset <bin offset>]\n"
          "  Calcuate the max/min coordinates (X,Y,Z) of atoms in <mask>.\n"
          "    [<mask>]: Atoms to calculate boundaries for.\n"
          "    [out <filename>]: Write boundaries to <filename>.\n"
          "    [dx <dx>] [dy <dy>] [dz <dz>]: Grid spacing in X/Y/Z directions.\n"
          "      If specified a grid will be created after processing is complete.\n"
          "    [name <gridname>]: Name of grid data set (if 'dx <dx>' etc specified).\n"
          "    [offset <bin offset>]: Number of bins to add in each direction to grid.\n");
}

// Action_Bounds::Init()
Action::RetType Action_Bounds::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  outfilename_ = actionArgs.GetStringKey("out");
  dxyz_[0] = actionArgs.getKeyDouble("dx", -1.0);
  dxyz_[1] = actionArgs.getKeyDouble("dy", -1.0);
  dxyz_[2] = actionArgs.getKeyDouble("dz", -1.0);
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  std::string dsname = actionArgs.GetStringKey("name");
  offset_ = actionArgs.getKeyInt("offset", 1);
  // Set default spacings if necessary. Y uses X, Z uses Y.
  if (dxyz_[0] > -1.0) {
    if (dsname.empty()) {
      mprinterr("Error: Grid name must be specified if spacing specified.\n");
      return Action::ERR;
    } 
    if (dxyz_[1] < 0.0) dxyz_[1] = dxyz_[0];
    if (dxyz_[2] < 0.0) dxyz_[2] = dxyz_[1];
    grid_ = DSL->AddSet(DataSet::GRID_FLT, dsname, "Bounds");
    if (grid_ == 0) return Action::ERR;
  }

  min_[0] = DBL_MAX;
  min_[1] = min_[0];
  min_[2] = min_[0];
  max_[0] = -DBL_MAX;
  max_[1] = max_[0];
  max_[2] = max_[0];

  mprintf("    BOUNDS: Calculating bounds for atoms in mask [%s]\n", mask_.MaskString());
  if (!outfilename_.empty())
    mprintf("\tOutput to file %s\n", outfilename_.c_str());
  if (grid_ != 0) {
    mprintf("\tGrid %s will be created after processing using\n"
            "\t  spacings dX= %g  dY= %g  dZ= %g  offset= %i Bins.\n",
            grid_->Legend().c_str(), dxyz_[0], dxyz_[1], dxyz_[2], offset_);
  }
  return Action::OK;
}

// Action_Bounds::Setup()
Action::RetType Action_Bounds::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask( mask_ ) ) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprinterr("Error: bounds: No atoms selected in mask.\n");
    return Action::ERR;
  }
  return Action::OK;
}

// Action_Bounds::DoAction()
Action::RetType Action_Bounds::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
  {
    const double* xyz = currentFrame->XYZ( *atom );
    if (xyz[0] < min_[0]) min_[0] = xyz[0];
    if (xyz[0] > max_[0]) max_[0] = xyz[0];
    if (xyz[1] < min_[1]) min_[1] = xyz[1];
    if (xyz[1] > max_[1]) max_[1] = xyz[1];
    if (xyz[2] < min_[2]) min_[2] = xyz[2];
    if (xyz[2] > max_[2]) max_[2] = xyz[2];
  }
  return Action::OK;
}

void Action_Bounds::Print() {
  static const char cXYZ[3] = {'X', 'Y', 'Z'};
  CpptrajFile outfile;
  Vec3 center;
  size_t nxyz[3];

  if ( outfile.OpenEnsembleWrite( outfilename_, ensembleNum_ ) ) return;
  for (int i = 0; i < 3; i++) {
    outfile.Printf("%f < %c < %f", min_[i], cXYZ[i], max_[i]);
    if (dxyz_[i] > 0.0) {
      center[i] = (max_[i] + min_[i]) / 2.0;
      long int nbins = (long int)ceil( (max_[i] - min_[i]) / dxyz_[i] ) + (long int)offset_;
      nxyz[i] = (size_t)nbins;
      outfile.Printf("\tCenter= %f  Bins=%zu", center[i], nxyz[i]);
    }
    outfile.Printf("\n");
  }
  outfile.CloseFile();
  if (grid_ != 0) {
    DataSet_3D& grid3d = static_cast<DataSet_3D&>( *grid_ );
    if (grid3d.Allocate_N_C_D( nxyz[0], nxyz[1], nxyz[2], center, dxyz_ ))
      mprinterr("Errror: Could not allocate grid %s\n", grid_->Legend().c_str());
  }
}
