#ifndef INC_DATASET_GRIDFLT_H
#define INC_DATASET_GRIDFLT_H
#include "DataSet_3D.h"
#include "Grid.h"
/// Single-precision three-dimensional grid.
class DataSet_GridFlt : public DataSet_3D {
  public:
    DataSet_GridFlt() : DataSet_3D(GRID_FLT, 12, 4) {}
    float& operator[](size_t idx)              { return grid_[idx];          }
    static DataSet* Alloc()       { return (DataSet*)new DataSet_GridFlt();  }
    Grid<float> const& InternalGrid()    const { return grid_; }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return grid_.size();        }
    int Sync()                                 { return 1;                   }
    void Info()                          const { return;                     }
    // ----- DataSet_3D functions ----------------
    int Allocate3D(size_t x,size_t y,size_t z) { return grid_.resize(x,y,z); }
    void Write3D(CpptrajFile&,int,int,int) const;
    double GetElement(int x,int y,int z) const { return (double)grid_.element(x,y,z); }
    void SetElement(int x,int y,int z,float v) { grid_.setGrid(x,y,z,v);     }
    double operator[](size_t idx)        const { return (double)grid_[idx];  }
    size_t NX() const { return grid_.NX(); }
    size_t NY() const { return grid_.NY(); }
    size_t NZ() const { return grid_.NZ(); }
    // -------------------------------------------
    /// Type definition of iterator over grid elements.
    typedef Grid<float>::iterator iterator;
    iterator begin() { return grid_.begin(); }
    iterator end()   { return grid_.end();   }
    /// Increment grid bin corresponding to point by given value.
    inline long int Increment(Vec3 const&, float);
    inline long int Increment(const double*, float);
    /// Increment grid bin by given value.
    inline long int Increment(int,int,int,float);
    /// \return grid value at specified bin.
    float GridVal(int x,int y,int z) const { return grid_.element(x,y,z); }
    /// \return grid index
    long int CalcIndex(int i, int j, int k) const {
      return grid_.CalcIndex(i,j,k);
    }
  private:
    Grid<float> grid_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// DataSet_GridFlt::Increment()
long int DataSet_GridFlt::Increment(Vec3 const& xyz, float f) {
  int i,j,k;
  if (CalcBins(xyz[0],xyz[1],xyz[2],i,j,k))
    return grid_.incrementBy(i,j,k,f);
  return -1L; 
}
// DataSet_GridFlt::Increment()
long int DataSet_GridFlt::Increment(const double* xyz, float f) {
  int i,j,k;
  if (CalcBins(xyz[0],xyz[1],xyz[2],i,j,k))
    return grid_.incrementBy(i,j,k,f);
  return -1L;
}
// DataSet_GridFlt::Increment()
long int DataSet_GridFlt::Increment(int i, int j, int k, float f) {
  return grid_.incrementBy(i,j,k,f);
}
#endif
