#ifndef INC_DATASET_COORDS_CRD_H
#define INC_DATASET_COORDS_CRD_H
#include "DataSet_Coords.h"
class DataSet_Coords_CRD : public DataSet_Coords {
  public:
    DataSet_Coords_CRD() : DataSet_Coords(COORDS) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_Coords_CRD();    }
    // ----- DataSet functions -------------------
    size_t Size() const { return coords_.size(); }
    int Sync()          { return 1;              }
    void Info() const;
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* )              { return;     }
    double Dval(size_t)                    const { return 0.0; }
    double Xcrd(size_t idx)  const { return Dim(0).Coord(idx); }
    void WriteBuffer(CpptrajFile&, size_t) const { return;     }
    // -------------------------------------------
    /// Add a frame.
    inline void AddFrame(Frame const& fIn) { 
      coords_.push_back( fIn.ConvertToCRD(numBoxCrd_, hasVel_) ); 
    }
    /// Get a frame at position.
    inline void GetFrame(int idx, Frame& fIn) { 
      fIn.SetFromCRD( coords_[idx], numCrd_, numBoxCrd_, hasVel_ ); 
    }
    /// Get a frame at position corresponding to mask.
    inline void GetFrame(int idx, Frame& fIn, AtomMask const& mIn) {
      fIn.SetFromCRD( coords_[idx], mIn, numCrd_, numBoxCrd_, hasVel_ );
    }
    /// Set CRD at position with frame.
    inline void SetCRD(int idx, Frame const& fIn) {
      coords_[idx] = fIn.ConvertToCRD(numBoxCrd_, hasVel_);
    }
  private:
    typedef std::vector<Frame::CRDtype> CRDarray;
    CRDarray coords_;                  ///< Array of coordinate frames.
};
#endif
