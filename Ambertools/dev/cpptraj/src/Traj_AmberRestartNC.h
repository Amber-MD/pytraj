#ifndef INC_TRAJ_AMBERRESTARTNC_H
#define INC_TRAJ_AMBERRESTARTNC_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
#include "NetcdfFile.h"
// Class: Traj_AmberRestartNC
/// TrajectoryIO class for reading and writing Amber Netcdf Restarts
class Traj_AmberRestartNC : public TrajectoryIO, private NetcdfFile {
  public:
    Traj_AmberRestartNC();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_AmberRestartNC(); }
    static void ReadHelp();
    static void WriteHelp();
    ~Traj_AmberRestartNC();
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    int processWriteArgs(ArgList&);
    int processReadArgs(ArgList&);
    void Info();
  private:
    double restartTime_;
    bool singleWrite_;
    bool useVelAsCoords_;
    bool outputTemp_;
    bool outputVel_;
    bool outputTime_;
    bool readAccess_;
    double time0_;
    double dt_;
    FileName filename_;

    int readVelocity(int, Frame&) { return 1; }
};
#endif
#endif  
