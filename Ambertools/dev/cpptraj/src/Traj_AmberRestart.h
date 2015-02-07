#ifndef INC_TRAJ_AMBERRESTART_H
#define INC_TRAJ_AMBERRESTART_H
#include "TrajectoryIO.h"
#include "BufferedFrame.h"
// Class: Traj_AmberRestart.h
/// Reads and writes formatted (ASCII text) amber
class Traj_AmberRestart : public TrajectoryIO {
  public:
    Traj_AmberRestart();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_AmberRestart(); }
    static void WriteHelp();
    static void ReadHelp();
  private:
    // Inherited functions
    int processReadArgs(ArgList&);
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int readVelocity(int, Frame&);
    int writeFrame(int,Frame const&);
    int processWriteArgs(ArgList&);
    void Info();

    int getBoxAngles(std::string const&, Box&);

    int restartAtoms_;     ///< Number of atoms in restart file
    int natom3_;           ///< Number of coords
    int numBoxCoords_;     ///< Number of box coords (3 or 6)
    size_t coordSize_;     ///< Size of coords in bytes, for reading past coords.
    double restartTime_;   ///< Time in restart file, read in
    double restartTemp_;   ///< (Optional) replica temperature, read in.
    double time0_;         ///< For writes, restart time offset
    double dt_;            ///< For writes, restart timestep (scaling)
    bool singleWrite_;     ///< If false, frame # will be appended to output filename
    bool readAccess_;      ///< If true, presence/absence of velocity info is known
    bool useVelAsCoords_;  ///< If true read velocities in as coordinates.
    bool outputTemp_;
    bool outputVel_;
    bool outputTime_;
    BufferedFrame file_;
};
#endif
