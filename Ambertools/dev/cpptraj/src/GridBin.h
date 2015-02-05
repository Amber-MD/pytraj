#ifndef INC_GRIDBIN_H
#define INC_GRIDBIN_H
#include "Matrix_3x3.h"
/// Class used to perform binning on/get voxel coords of 3D grids.
class GridBin {
  public:
    GridBin() : OXYZ_(0.0) {}
    virtual ~GridBin() {}
    /// Given coordinates, set corresponding bin indices; check bounds.
    virtual bool CalcBins(double, double, double, int&, int&, int&) const = 0;
    /// Given coordinates, set corresponding bin indices; no bounds check.
    virtual void BinIndices(double, double, double, int&, int&, int&) const = 0;
    /// \return coordinates of bin for given indices.
    virtual Vec3 BinCorner(int, int, int) const = 0;
    /// \return coordinates of bin center for given indices.
    virtual Vec3 BinCenter(int, int, int) const = 0;
    /// \return unit cell matrix. // TODO: Make const&?
    virtual Matrix_3x3 Ucell() const = 0;
    /// \return true if GridBin type is orthogonal.
    virtual bool IsOrthoGrid() const = 0;
    /// \return Voxel volume.
    virtual double VoxelVolume() const = 0;
    /// \return Grid origin.
    Vec3 const& GridOrigin() const { return OXYZ_; }
  protected:
    Vec3 OXYZ_; ///< Grid origin.
};
// -----------------------------------------------------------------------------
/// Orthogonal grid
class GridBin_Ortho : public GridBin {
  public:
    GridBin_Ortho() : dx_(-1.0), dy_(-1.0), dz_(-1.0), mx_(0),  my_(0), mz_(0) {}
    bool CalcBins(double x, double y, double z,
                  int& i, int& j, int& k) const
    {
      if (x >= OXYZ_[0] && x < mx_) { // X
        if (y >= OXYZ_[1] && y < my_) { // Y
          if (z >= OXYZ_[2] && z < mz_) { // Z
            i = (int)((x-OXYZ_[0]) / dx_);
            j = (int)((y-OXYZ_[1]) / dy_);
            k = (int)((z-OXYZ_[2]) / dz_);
            return true;
          }
        }
      }
      return false;
    }
    void BinIndices(double x, double y, double z, int& i, int& j, int& k) const {
      i = (int)((x-OXYZ_[0]) / dx_);
      j = (int)((y-OXYZ_[1]) / dy_);
      k = (int)((z-OXYZ_[2]) / dz_);
    }
    Vec3 BinCorner(int i, int j, int k) const {
      return Vec3((double)i*dx_+OXYZ_[0],
                  (double)j*dy_+OXYZ_[1],
                  (double)k*dz_+OXYZ_[2]);
    }
    Vec3 BinCenter(int i, int j, int k) const {
      return Vec3((double)i*dx_+OXYZ_[0]+0.5*dx_,
                  (double)j*dy_+OXYZ_[1]+0.5*dy_,
                  (double)k*dz_+OXYZ_[2]+0.5*dz_);
    }
    Matrix_3x3 Ucell() const { return Matrix_3x3(mx_-OXYZ_[0], my_-OXYZ_[1], mz_-OXYZ_[2]); }
    bool IsOrthoGrid() const { return true; }
    double VoxelVolume() const { return dx_ * dy_ * dz_; }
    /// Setup with given origin, spacing; calculate maximum.
    void Setup_O_D(size_t nx, size_t ny, size_t nz,
                   Vec3 const& oxyzIn, Vec3 const& dxyz)
    {
      OXYZ_ = oxyzIn;
      dx_ = dxyz[0]; dy_ = dxyz[1]; dz_ = dxyz[2];
      mx_ = OXYZ_[0] + ((double)nx * dx_);
      my_ = OXYZ_[1] + ((double)ny * dy_);
      mz_ = OXYZ_[2] + ((double)nz * dz_);
    }
    double DX() const { return dx_; }
    double DY() const { return dy_; }
    double DZ() const { return dz_; }
  private:
    double dx_, dy_, dz_; ///< Grid spacing
    double mx_, my_, mz_; ///< Grid max
};
// -----------------------------------------------------------------------------
/// Non-orthogonal grid.
class GridBin_Nonortho : public GridBin {
  public:
    GridBin_Nonortho() {}
    bool CalcBins(double x, double y, double z,
                  int& i, int& j, int& k) const
    {
      Vec3 frac = recip_ * Vec3(x - OXYZ_[0], y - OXYZ_[1], z - OXYZ_[2]);
      if (frac[0] >= 0.0 && frac[0] < 1.0) {
        if (frac[1] >= 0.0 && frac[1] < 1.0) {
          if (frac[2] >= 0.0 && frac[2] < 1.0) {
            i = (int)(frac[0] * nx_);
            j = (int)(frac[1] * ny_);
            k = (int)(frac[2] * nz_);
            return true;
          }
        }
      }
      return false;
    }
    void BinIndices(double x, double y, double z, int& i, int& j, int& k) const {
      Vec3 frac = recip_ * Vec3(x - OXYZ_[0], y - OXYZ_[1], z - OXYZ_[2]);
      i = (int)(frac[0] * nx_);
      j = (int)(frac[1] * ny_);
      k = (int)(frac[2] * nz_);
    }
    Vec3 BinCorner(int i, int j, int k) const {
      Vec3 frac( (double)i / nx_, (double)j / ny_, (double)k / nz_ );
      return ucell_.TransposeMult( frac );
    }
    Vec3 BinCenter(int i, int j, int k) const {
      Vec3 frac_half((1.0 + 2.0 * (double)i) / (2.0 * nx_),  //(0.5 * (1.0 / nx_)) + ((double)i / nx_),
                     (1.0 + 2.0 * (double)j) / (2.0 * ny_), 
                     (1.0 + 2.0 * (double)k) / (2.0 * nz_));
      return ucell_.TransposeMult( frac_half );
    }
    Matrix_3x3 Ucell() const { return ucell_; }
    bool IsOrthoGrid() const { return false; }
    double VoxelVolume() const { return voxelvolume_; }
    /// Setup with given bins, origin and box coordinates.
    void Setup_O_Box(size_t nxIn, size_t nyIn, size_t nzIn,
                     Vec3 const& oxyzIn, Box const& boxIn) {
      nx_ = (double)nxIn; ny_ = (double)nyIn; nz_ = (double)nzIn;
      OXYZ_ = oxyzIn;
      // Get unit cell and fractional cell vectors (recip).
      double vol = boxIn.ToRecip( ucell_, recip_ );
      voxelvolume_ = vol / (nx_ * ny_ * nz_);
    }
  private:
    double nx_, ny_, nz_; ///< Number of bins in double precision.
    double voxelvolume_;  ///< Voxel volume.
    Matrix_3x3 ucell_;    ///< Unit cell axes coordinates.
    Matrix_3x3 recip_;    ///< Fractional axes coordinates.
};
#endif
