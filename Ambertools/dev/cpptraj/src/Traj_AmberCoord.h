#ifndef INC_TRAJ_AMBERCOORD_H
#define INC_TRAJ_AMBERCOORD_H
#include "TrajectoryIO.h"
#include "BufferedFrame.h"
// Class: Traj_AmberCoord
/// Reads and writes formatted (ASCII text) amber trajectories. 
class Traj_AmberCoord: public TrajectoryIO {
  public:
    Traj_AmberCoord();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_AmberCoord(); }
    static void WriteHelp();
  private:
    int natom3_;          ///< Number of coords (# atoms X 3)
    // TODO: Replace hasREMD with REMD_HEADER_SIZE/HasT()
    size_t headerSize_;   ///< REMD header size if present
    size_t tStart_;       ///< Start position for T-value in header
    size_t tEnd_;         ///< End position for T-value in header
    int numBoxCoords_;    ///< Number of box coords, 3 (ortho or truncoct) or 6 (triclinic)
    const char* outfmt_;  ///< Format string for writing coordinates
    bool highPrecision_;  ///< If true output format will be 8.6 instead of 8.3
    bool outputTemp_;     ///< If true temperature will be out to REMD line.
    BufferedFrame file_;  ///< Buffer reads/writes
    double boxAngle_[3];  ///< Hold default box angles in case traj has only box lengths

    static const size_t REMD_HEADER_SIZE;
    static const size_t RXSGLD_HEADER_SIZE;

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&);
    int readVelocity(int, Frame&);
    int processReadArgs(ArgList&)  { return 0; }
};
#endif
