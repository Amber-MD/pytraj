#ifndef INC_TRAJ_AMBERNETCDF_H
#define INC_TRAJ_AMBERNETCDF_H
#ifdef BINTRAJ
#include "TrajectoryIO.h"
#include "NetcdfFile.h"
// Class: Traj_AmberNetcdf
/// Reads and writes Amber Netcdf format trajectories. 
class Traj_AmberNetcdf : public TrajectoryIO, private NetcdfFile {
  public:
    Traj_AmberNetcdf();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_AmberNetcdf(); }
    ~Traj_AmberNetcdf();
    static void ReadHelp();
    static void WriteHelp();
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int readVelocity(int, Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&);
    int processReadArgs(ArgList&);
    // Reservoir functions
    inline int createReservoir(bool,double,int);
    int writeReservoir(int, Frame&, double, int);
  private:
    float *Coord_;
    FileName filename_;
    int eptotVID_;
    int binsVID_;
    bool useVelAsCoords_;
    bool readAccess_;
    bool outputTemp_;
    bool outputVel_;
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
int Traj_AmberNetcdf::createReservoir(bool hasBins, double reservoirT, int iseed) {
  return NC_createReservoir(hasBins, reservoirT, iseed, eptotVID_, binsVID_);
}
#endif
#endif
