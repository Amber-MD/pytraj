#ifndef INC_TRAJ_MOL2FILE_H
#define INC_TRAJ_MOL2FILE_H
#include "TrajectoryIO.h"
#include "Mol2File.h"
/// TrajectoryIO class for reading coordinates from Mol2 files.
class Traj_Mol2File : public TrajectoryIO {
  public:
    /// Indicate how the mol2 file should be written.
    /** - SINGLE: Writing only a single frame
      * - MOL: Multiple frames written to the same file separated with
      *        a @<TRIPOS>MOLECULE section.
      * - MULTI: Each frame written to a different file with name filename.frame
      */
    enum MOL2WRITEMODE { NONE = 0, SINGLE, MOL, MULTI };

    Traj_Mol2File();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_Mol2File(); }
    static void WriteHelp();
  private:
    MOL2WRITEMODE mol2WriteMode_;
    Topology* mol2Top_;
    int currentSet_;
    bool hasCharges_;
    Mol2File file_; 

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
