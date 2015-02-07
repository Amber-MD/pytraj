#ifndef INC_TRAJ_SDF_H
#define INC_TRAJ_SDF_H
#include "TrajectoryIO.h"
#include "SDFfile.h"
/// Read coordinates from SDF file. Limited to 1 frame currently.
class Traj_SDF : public TrajectoryIO {
  public:
    Traj_SDF() {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_SDF(); }
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int readVelocity(int, Frame&) { return 1; }
    int processWriteArgs(ArgList&) { return 0; }
    int processReadArgs(ArgList&) { return 0; }

    SDFfile file_;
};
#endif
