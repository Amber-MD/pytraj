#ifndef INC_TRAJ_CHARMMDCD_H
#define INC_TRAJ_CHARMMDCD_H
#include "TrajectoryIO.h"
// Class: Traj_CharmmDcd
/// TrajectoryIO class for reading coordinates from charmm dcd files.
class Traj_CharmmDcd : public TrajectoryIO {
  public :
    Traj_CharmmDcd();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_CharmmDcd(); }
    static void WriteHelp();
    ~Traj_CharmmDcd();
  private:
    int dcdatom_;            ///< Number of atoms in DCD file.
    int dcdframes_;          ///< Number of frames in DCD file.
    bool isBigEndian_;       ///< True if file is Big endian
    bool is64bit_;           ///< True if file is 64 bit
    unsigned int blockSize_; ///< Size of block bytes: 32 bit = 4, 64 bit = 8
    size_t dcd_dim_;         ///< Number of dimensions in DCD file.
    size_t boxBytes_;        ///< Number of bytes used by box coords if present.
    size_t frame1Bytes_;     ///< Number of bytes used by first frame.
    size_t frameNBytes_;     ///< Number of bytes used by other frames (==frame1 if namnf==0).
    size_t headerBytes_;     ///< Size of DCD header in bytes.
    size_t coordinate_size_; ///< Size of X|Y|Z coord frame in bytes.
    int nfixedat_;           ///< Number of fixed atoms
    int nfreeat_;            ///< Number of free atoms
    int* freeat_;            ///< Free atom indices
    float* xcoord_;          ///< Master coord array, start of X coords
    float* ycoord_;          ///< Pointer to start of Y coords in master coord array
    float* zcoord_;          ///< Pointer to start of Z coords in master coord array
    CpptrajFile file_;       ///< Input/Output file

    union headerbyte { unsigned char c[80]; int i[20]; float f[20]; };
    int ReadBlock(int);
    int WriteBlock(int);
    void AllocateCoords();
    int readDcdHeader();
    int ReadBox(double*);
    int writeDcdHeader();

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

    int readVelocity(int, Frame&) { return 1; }
    int processReadArgs(ArgList&) { return 0; }
};
#endif
