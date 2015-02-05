#ifndef INC_REFERENCEFRAME_H
#define INC_REFERENCEFRAME_H
#include "DataSet_Coords_REF.h"
/// Wrapper around DataSet_Coords_REF DataSet.
/** Intended as a non-modifiable holder for the DataSet_Coords_REF DataSet
  * that also can hold an error status.
  */
class ReferenceFrame {
  public:
    ReferenceFrame() : ref_(0), err_(0) {}
    ReferenceFrame(int err) : ref_(0), err_(err) {}
    ReferenceFrame(DataSet_Coords_REF* ds) : ref_(ds), err_(0) {}
    ReferenceFrame(const ReferenceFrame& rhs) : ref_(rhs.ref_), err_(rhs.err_) {}
    Frame const& Coord()        const { return ref_->RefFrame();  }
    Topology const& Parm()      const { return ref_->Top();       }
    bool error()                const { return err_ != 0;         }
    bool empty()                const { return ref_ == 0;         }
    FileName const& FrameName() const { return ref_->FrameName(); }
    DataSet_Coords_REF* RefPtr()const { return ref_;              }
  private:
    DataSet_Coords_REF* ref_; ///< Reference coords DataSet from e.g. DataSetList.
    int err_;                 ///< Error status.
};
#endif
