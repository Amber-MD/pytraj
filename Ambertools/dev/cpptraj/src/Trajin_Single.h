#ifndef TRAJIN_SINGLE_H
#define TRAJIN_SINGLE_H
#include "Trajin.h"
/// Class for reading in single trajectories.
class Trajin_Single : public Trajin {
  public:
    Trajin_Single();
    ~Trajin_Single();
    /// Set up trajectory for reading, optionally checking box info.
    int SetupTrajRead(std::string const&, ArgList&, Topology*, bool);
    /// Set up trajectory for reading, check box info.
    int SetupTrajRead(std::string const&, ArgList&, Topology*);
    int BeginTraj(bool);
    void EndTraj();
    int ReadTrajFrame(int, Frame&);
    void PrintInfo(int) const;
    CoordinateInfo const& TrajCoordInfo() const { return cInfo_; }
    // -------------------------------------------
    std::string Title() {
      if (trajio_==0) return std::string("");
      else return trajio_->Title();
    }
  private:
    TrajectoryIO* trajio_; ///< Hold class that will interface with traj format.
    TrajectoryIO* velio_;  ///< Hold class that will interface with opt. mdvel file.
    CoordinateInfo cInfo_; ///< Hold coordinate metadata.
    bool trajIsOpen_;      ///< True if trajectory is open. 
};
#endif
