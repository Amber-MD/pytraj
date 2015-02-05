#ifndef INC_TRAJ_PDBFILE_H
#define INC_TRAJ_PDBFILE_H
#include "TrajectoryIO.h"
#include "PDBfile.h"
// Class: Traj_PDBfile
/// TrajectoryIO class for reading coordinates from PDB files.
class Traj_PDBfile: public TrajectoryIO {
  public:
    /** PDBWRITEMODE: Indicate how the pdb should be written.
      *  SINGLE: Writing only a single frame.
      *  MODEL: Multiple frames written to the same file separated with 
      *         the MODEL keyword
      *  MULTI: Each frame written to a different file with name filename.frame
      */
    enum PDBWRITEMODE {NONE = 0, SINGLE, MODEL, MULTI};

    Traj_PDBfile();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_PDBfile(); }
    static void WriteHelp();
  private:
    int pdbAtom_;
    int currentSet_;
    int ter_num_; ///< Amount to increment atom number for TER
    PDBWRITEMODE pdbWriteMode_;
    bool dumpq_;   ///< If true print charges in Occupancy column
    bool dumpr_;   ///< If true print radii in B-factor column.
    bool pdbres_;  ///< If true convert Amber res names to PDBV3 style.
    bool pdbatom_; ///< If true convert Amber atom names to PDBV3 style.
    bool write_cryst1_; ///< If false write CRYST1 for first frame.
    std::string space_group_;
    Topology *pdbTop_;
    PDBfile file_;

    std::vector<char> chainID_;      ///< Hold chainID for each atom.
    std::vector<NameType> resNames_; ///< Hold residue names.
    char chainchar_;

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(std::string const&, Topology*);
    int setupTrajout(std::string const&, Topology*, int, bool);
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
