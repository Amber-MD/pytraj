#ifndef INC_DATASET_3D_H
#define INC_DATASET_3D_H
#include "DataSet.h"
#include "CpptrajFile.h"
#include "Box.h"
#include "GridBin.h"
/// Interface for 3D DataSets.
// FIXME: Use DataSet Dims?
class DataSet_3D : public DataSet {
  public:
    DataSet_3D() : gridBin_(0) {}
    virtual ~DataSet_3D(); // Virtual since this class is inherited.
    DataSet_3D(DataSet::DataType tIn, int wIn, int pIn) :
      DataSet(tIn, wIn, pIn, 3), gridBin_(0) {}
    /// Write 3D data to file.
    virtual void Write3D(CpptrajFile&,int,int,int) const = 0;
    /// \return Data from grid at x/y/z point.
    virtual double GetElement(int, int, int) const = 0;
    /// \return Data from grid.
    virtual double operator[](size_t) const = 0;
    /// \return size of X dimension.
    virtual size_t NX() const = 0;
    /// \return size of Y dimension.
    virtual size_t NY() const = 0;
    /// \return size of Z dimension.
    virtual size_t NZ() const = 0;
    // -------------------------------------------
    // TODO: Remove this. Only needed by DataSet_1D.h
    void Add(size_t,const void*) { }
    /// Set up grid from dims, origin, and spacing.
    int Allocate_N_O_D(size_t,size_t,size_t,Vec3 const&,Vec3 const&);
    /// Set up grid from dims, center, and spacing.
    int Allocate_N_C_D(size_t,size_t,size_t,Vec3 const&,Vec3 const&);
    /// Set up grid from sizes, center, and spacing.
    int Allocate_X_C_D(Vec3 const&,Vec3 const&,Vec3 const&);
    /// Set up grid from dims, origin, and box.
    int Allocate_N_O_Box(size_t,size_t,size_t, Vec3 const&, Box const&);
    /// Print grid info.
    void GridInfo() const;
    /// Convert X, Y, and Z coords to indices. Check bounds.
    bool CalcBins(double x,double y,double z,int& i,int& j,int& k) const { 
      return gridBin_->CalcBins(x, y, z, i, j, k);
    }
    /// Convert X, Y, and Z coords to indices. No bounds check.
    void BinIndices(double x,double y,double z,int& i,int& j,int& k) const {
      return gridBin_->BinIndices(x, y, z, i, j, k);
    }
    /// \return coordinates of specified voxel corner.
    Vec3 BinCorner(int i,int j,int k) const {
      return gridBin_->BinCorner(i, j, k);
    }
    /// \return coordinates of specified voxel center.
    Vec3 BinCenter(int i,int j,int k) const {
      return gridBin_->BinCenter(i, j, k);
    }
    /// \return coordinates of grid origin.
    Vec3 const& GridOrigin() const {
      return gridBin_->GridOrigin();
    }
    /// \return unit cell matrix.
    Matrix_3x3 Ucell() const { return gridBin_->Ucell(); }
    /// \return voxel volume.
    double VoxelVolume() const { return gridBin_->VoxelVolume(); }
  private:
    /// Check if grid dimension is even; if not, increment it by 1.
    static void CheckEven(size_t&, char);
    /// Set up grid for given # x, y, and z points.
    // TODO: Make public if grids will be used for other than binning.
    virtual int Allocate3D(size_t, size_t, size_t) = 0;

    GridBin* gridBin_; ///< Used to calculate bins/coords depending on grid type.
};
#endif
