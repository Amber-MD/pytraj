#ifndef INC_DATASET_COORDS_H
#define INC_DATASET_COORDS_H
#include "Topology.h"
#include "DataSet_1D.h"
// NOTE: Should this be a 4D DataSet?
/// Interface to COORDS data sets
class DataSet_Coords : public DataSet_1D {
  public:
    DataSet_Coords() : numCrd_(0), numBoxCrd_(0), hasVel_(false) {}
    DataSet_Coords(DataSet::DataType t) : 
      DataSet_1D(t, 8, 3), numCrd_(0), numBoxCrd_(0), hasVel_(false) {}
    virtual ~DataSet_Coords() {}
    /// Allocate a Frame that can be used to store COORDS 
    Frame AllocateFrame() const;
    /// Add given Frame to this COORDS
    virtual void AddFrame(Frame const&) = 0;
    /// Set COORDS at specified position with Frame
    virtual void SetCRD(int, Frame const&) = 0;
    /// Set given Frame with COORDS at specified position
    virtual void GetFrame(int, Frame&) = 0;
    /// Set given Frame with COORDS at specified position according to mask
    virtual void GetFrame(int, Frame&, AtomMask const&) = 0;
    /// Set main topology that will be associated with frames to/from this COORDS
    void SetTopology(Topology const&);
    inline Topology const& Top() const { return top_; }
  protected:
    // TODO: Make unsigned
    int numCrd_;    ///< Number of coordinates
    int numBoxCrd_; ///< Number of box coords (0 or 6).
    bool hasVel_;   ///< Coordinates contain velocities.
    Topology top_;  ///< Topology corresponding to coordinates.
};
#endif
