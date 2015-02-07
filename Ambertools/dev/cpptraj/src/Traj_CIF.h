#ifndef INC_TRAJ_CIF_H
#define INC_TRAJ_CIF_H
#include "TrajectoryIO.h"
#include "CIFfile.h"
// Class: Traj_CIF
/// TrajecttoryIO class for reading coordinates from CIF files.
class Traj_CIF : public TrajectoryIO {
  public:
    Traj_CIF() : Natoms_(0), Nmodels_(0), Cartn_x_col_(0),
                 Cartn_y_col_(0), Cartn_z_col_(0) {}
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_CIF(); }
  private:
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int openTrajin();
    int readFrame(int,Frame&);
    void Info();
    void closeTraj() {}
    int processWriteArgs(ArgList&) { return 0; }
    int readVelocity(int, Frame&)  { return 1; }
    int processReadArgs(ArgList&)  { return 0; }
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool) { return 1; }
    int writeFrame(int,Frame const&)                           { return 1; } 

    CIFfile file_;
    Box boxInfo_;
    int Natoms_;
    int Nmodels_;
    int Cartn_x_col_;
    int Cartn_y_col_;
    int Cartn_z_col_;
};
#endif
