#ifndef INC_DATASET_COORDS_H
#define INC_DATASET_COORDS_H
#include "DataSet.h"
#include "Topology.h"
/// Interface to COORDS data sets
class DataSet_Coords : public DataSet {
  public:
    DataSet_Coords() : numCrd_(0), numBoxCrd_(0), hasVel_(false) {}
    DataSet_Coords(DataSet::DataType t) : 
      DataSet(t, 8, 3, 4), numCrd_(0), numBoxCrd_(0), hasVel_(false) {}
    virtual ~DataSet_Coords() {}
    /// Allocate memory for a certain number of frames.
    virtual int AllocateCoords(size_t) = 0;
    /// Add given Frame to this COORDS
    virtual void AddFrame(Frame const&) = 0;
    /// Set COORDS at specified position with Frame
    virtual void SetCRD(int, Frame const&) = 0;
    /// Set given Frame with COORDS at specified position
    virtual void GetFrame(int, Frame&) = 0;
    /// Set given Frame with COORDS at specified position according to mask
    virtual void GetFrame(int, Frame&, AtomMask const&) = 0;
    // -------------------------------------------
    /// Allocate a Frame that can be used to store COORDS 
    Frame AllocateFrame() const;
    /// Set main topology that will be associated with frames to/from this COORDS
    void SetTopology(Topology const&);
    /// \return topology associated with these COORDS.
    inline Topology const& Top() const { return top_; }
  protected:
    // TODO: Make unsigned
    int numCrd_;    ///< Number of coordinates
    int numBoxCrd_; ///< Number of box coords (0 or 6).
    bool hasVel_;   ///< Coordinates contain velocities.
    Topology top_;  ///< Topology corresponding to coordinates.
};
#endif
