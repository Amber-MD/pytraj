#ifndef INC_TRAJ_BINPOS
#define INC_TRAJ_BINPOS
#include "TrajectoryIO.h"
// Class: Traj_Binpos
/// Read and write SCRIPPS BINPOS format
class Traj_Binpos : public TrajectoryIO {
  public:
    Traj_Binpos();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_Binpos(); }
    ~Traj_Binpos();
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
    int processWriteArgs(ArgList&) { return 0; }

    int bpnatom_;
    int bpnatom3_;
    size_t frameSize_;
    float* bpbuffer_;
    CpptrajFile file_;

    int readVelocity(int, Frame&) { return 1; }
    int processReadArgs(ArgList&) { return 0; }
};
#endif
