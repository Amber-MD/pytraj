#ifndef INC_COORDINATEINFO_H
#define INC_COORDINATEINFO_H
#include "ReplicaDimArray.h"
#include "Box.h"
/// All metadata associated with a Frame.
class CoordinateInfo {
  public:
    /// CONSTRUCTOR
    CoordinateInfo() : hasVel_(false), hasTemp_(false), hasTime_(false), hasFrc_(false) {}
    /// CONSTRUCTOR - box, velocity, temperature, time
    CoordinateInfo(Box const& b, bool v, bool t, bool m) :
      box_(b), hasVel_(v), hasTemp_(t), hasTime_(m), hasFrc_(false) {}
    /// CONSTRUCTOR - all
    CoordinateInfo(ReplicaDimArray const& r, Box const& b, bool v, bool t, bool m, bool f) :
      remdDim_(r), box_(b), hasVel_(v), hasTemp_(t), hasTime_(m), hasFrc_(f) {}
    bool HasBox()              const { return box_.HasBox();            }
    const Box& TrajBox()       const { return box_;                     }
    bool HasVel()              const { return hasVel_;                  }
    bool HasTemp()             const { return hasTemp_;                 }
    bool HasTime()             const { return hasTime_;                 }
    bool HasForce()            const { return hasFrc_;                  }
    bool HasReplicaDims()      const { return (remdDim_.Ndims() != 0);  }
    ReplicaDimArray const& ReplicaDimensions() const { return remdDim_; }
    void SetTime(bool m)        { hasTime_ = m; }
    void SetTemperature(bool t) { hasTemp_ = t; }
    void SetVelocity(bool v)    { hasVel_ = v;  }
    void SetBox(Box const& b)   { box_ = b;     }
  private:
    ReplicaDimArray remdDim_; ///< Hold info on any replica dimensions.
    Box box_;                 ///< Hold box information.
    bool hasVel_;             ///< True if coords have associated velocities.
    bool hasTemp_;            ///< True if coords include temp info.
    bool hasTime_;            ///< True if coords include time info.
    bool hasFrc_;             ///< True if coords have associated forces.
};
#endif
