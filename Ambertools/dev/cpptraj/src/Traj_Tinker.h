#ifndef INC_TRAJ_TINKER_H
#define INC_TRAJ_TINKER_H
#include "TrajectoryIO.h"
#include "TinkerFile.h"
/// TrajectoryIO class for reading coordinates from Tinker XYZ/ARC files.
class Traj_Tinker : public TrajectoryIO {
  public:
    Traj_Tinker();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_Tinker(); }
  private:
    Topology* tinkerTop_;
    int currentSet_;
    TinkerFile file_; 

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&) { return 0; }
    int readVelocity(int, Frame&)  { return 1; }
    int processReadArgs(ArgList&)  { return 0; }
};
#endif
