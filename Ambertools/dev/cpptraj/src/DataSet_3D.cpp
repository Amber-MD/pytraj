#include "DataSet_3D.h"
#include "CpptrajStdio.h"

// DESTRUCTOR
DataSet_3D::~DataSet_3D() { if (gridBin_ != 0) delete gridBin_; }

// DataSet_3D::Allocate_N_O_Box()
int DataSet_3D::Allocate_N_O_Box(size_t nx, size_t ny, size_t nz, 
                                 Vec3 const& oxyz, Box const& boxIn)
{
  if (nx == 0 || ny == 0 || nz == 0) return 1;
  if (gridBin_ != 0) delete gridBin_;
  GridBin_Nonortho* gb = new GridBin_Nonortho();
  // Set origin and unit cell params.
  gb->Setup_O_Box(nx, ny, nz, oxyz, boxIn);
  gridBin_ = (GridBin*)gb;
  return Allocate3D(nx, ny, nz);
}

// DataSet_3D::Allocate_N_O_D()
int DataSet_3D::Allocate_N_O_D(size_t nx, size_t ny, size_t nz,
                               Vec3 const& oxyz, Vec3 const& dxyz)
{
  if (nx == 0 || ny == 0 || nz == 0) return 1;
  if (gridBin_ != 0) delete gridBin_;
  GridBin_Ortho* gb = new GridBin_Ortho();
  // Set origin and spacing, calculate maximum (for binning).
  gb->Setup_O_D(nx, ny, nz, oxyz, dxyz);
  gridBin_ = (GridBin*)gb;
  return Allocate3D(nx, ny, nz);
}

// Calc_Origin()
/** For even-spaced grids, origin is center - (N/2)*spacing.
  * For odd-spaced grids, origin is center - ((N-1/2)*spacing)+half_spacing
  */
static double Calc_Origin(int N, double D) {
  int odd = N % 2;
  int half = (N - odd) / 2;
  if (odd)
    return -(((double)half * D) + (D * 0.5));
  else
    return -((double)half * D);
}

// DataSet_3D::Allocate_N_C_D()
int DataSet_3D::Allocate_N_C_D(size_t nx, size_t ny, size_t nz,
                               Vec3 const& cxyz, Vec3 const& dxyz)
{
  // Calculate origin from center coordinates.
  Vec3 oxyz( cxyz[0] + Calc_Origin(nx, dxyz[0]),
             cxyz[1] + Calc_Origin(ny, dxyz[1]),
             cxyz[2] + Calc_Origin(nz, dxyz[2]) );
  return Allocate_N_O_D(nx,ny,nz,oxyz,dxyz);
}

// DataSet_3D::Allocate_X_C_D()
int DataSet_3D::Allocate_X_C_D(Vec3 const& sizes, Vec3 const& center, Vec3 const& dxyz)
{
  // Calculate bin counts
  size_t nx = (size_t)(sizes[0] / dxyz[0]);
  size_t ny = (size_t)(sizes[1] / dxyz[1]);
  size_t nz = (size_t)(sizes[2] / dxyz[2]);
  return Allocate_N_C_D( nx, ny, nz, center, dxyz );
}

void DataSet_3D::GridInfo() const {
  if (gridBin_ == 0) return;
  Vec3 const& oxyz = gridBin_->GridOrigin();
  mprintf("\t\t-=Grid Dims=- %8s %8s %8s\n", "X", "Y", "Z");
  mprintf("\t\t        Bins: %8zu %8zu %8zu\n", NX(), NY(), NZ());
  mprintf("\t\t      Origin: %8g %8g %8g\n", oxyz[0], oxyz[1], oxyz[2]);
  if (gridBin_->IsOrthoGrid()) {
    GridBin_Ortho const& gb = static_cast<GridBin_Ortho const&>( *gridBin_ );
    mprintf("\t\t     Spacing: %8g %8g %8g\n", gb.DX(), gb.DY(), gb.DZ());
    mprintf("\t\t      Center: %8g %8g %8g\n",
            oxyz[0] + (NX()/2)*gb.DX(),
            oxyz[1] + (NY()/2)*gb.DY(),
            oxyz[2] + (NZ()/2)*gb.DZ());
    //mprintf("\tGrid max    : %8.3f %8.3f %8.3f\n", grid.MX(), grid.MY(), grid.MZ());
  } else {
    Box box(gridBin_->Ucell());
    mprintf("\t\tBox: %s ABC={%g %g %g} abg={%g %g %g}\n", box.TypeName(),
            box[0], box[1], box[2], box[3], box[4], box[5]);
  }
}
