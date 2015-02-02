#ifndef INC_PARM_AMBER_H
#define INC_PARM_AMBER_H
#include "ParmIO.h"
class Parm_Amber : public ParmIO {
  public :
    Parm_Amber();
    ~Parm_Amber();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_Amber(); }
    static void WriteHelp();
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(std::string const&, Topology&);
    int processWriteArgs(ArgList&);
    int WriteParm(std::string const&, Topology const&);
    void SetDebug(int debugIn) { debug_ = debugIn; }
  private :
    /// Enumerated type for Fortran data type
    enum FortranType {
      UNKNOWN_FTYPE=0, FINT, FDOUBLE, FCHAR, FFLOAT
    };
    /// Enumerated type for Amber Parmtop Flags
    enum AmberParmFlagType {
      F_POINTERS = 0, F_NAMES,     F_CHARGE,    F_MASS,     F_RESNAMES,
      F_RESNUMS,      F_TYPES,     F_BONDSH,    F_BONDS,    F_SOLVENT_POINTER,
      F_ATOMSPERMOL,  F_PARMBOX,   F_ATYPEIDX,  F_NUMEX,    F_NB_INDEX,
      F_LJ_A,         F_LJ_B,      F_EXCLUDE,   F_RADII,    F_SCREEN,
      F_BONDRK,       F_BONDREQ,   F_ANGLETK,   F_ANGLETEQ, F_DIHPK,
      F_DIHPN,        F_DIHPHASE,  F_SCEE,      F_SCNB,     F_SOLTY,
      F_ANGLESH,      F_ANGLES,    F_DIHH,      F_DIH,      F_ASOL,
      F_BSOL,         F_HBCUT,     F_ITREE,     F_JOIN,     F_IROTAT,
      F_ATOMICNUM,    F_TITLE,     F_RADSET,    F_LES_NTYP, F_LES_TYPE,
      F_LES_FAC,      F_LES_CNUM,  F_LES_ID,    F_CAP_INFO, F_CAP_INFO2,
      F_IPOL,         F_POLAR,     F_CTITLE,    F_CHM_UBC,  F_CHM_UB,
      F_CHM_UBFC,     F_CHM_UBEQ,  F_CHM_NIMP,  F_CHM_IMP,  F_CHM_NIMPT,
      F_CHM_IMPFC,    F_CHM_IMPP,  F_LJ14A,     F_LJ14B,    F_CHM_CMAPC,
      F_CHM_CMAPR,    F_CHM_CMAPP, F_CHM_CMAPI, F_FF_TYPE,  F_PDB_RES,
      F_PDB_CHAIN,    F_PDB_ICODE, F_PDB_ALT
    };
    static const int AMBERPOINTERS;
    struct ParmFlag {
      const char* Flag; ///< %FLAG name in topology.
      const char* Fmt;  ///< Fortran format string for writing.
    };
    static const ParmFlag FLAGS[];

    // NOTE: Although amber topology files should only ever require 83 chars
    //       to read each line (80 chars + newline + CR (if dos) + NULL)
    //       increase the size to handle non-standard files.
    static const size_t BUF_SIZE = 256;
    char lineBuffer_[BUF_SIZE];
    int debug_;
    bool nochamber_; ///< For writes when true do not print chamber info
    enum ParmType { OLDPARM = 0, NEWPARM, CHAMBER };
    ParmType ptype_;
    std::string fformat_;
    FortranType ftype_;
    int fncols_;
    int fprecision_;
    int fwidth_;
    int error_count_;
    char *buffer_;
    size_t buffer_size_;
    size_t buffer_max_size_;
    CpptrajFile file_;

    int ReadAmberParm(Topology&);

    std::string GetLine();
    std::vector<int> GetInteger(int,int,int);
    std::vector<double> GetDouble(int,int,int);
    std::vector<NameType> GetName(int,int,int);
    std::vector<int> GetFlagInteger(AmberParmFlagType,int);
    std::vector<double> GetFlagDouble(AmberParmFlagType,int);
    std::vector<NameType> GetFlagName(AmberParmFlagType,int);
    bool SeekToFlag(AmberParmFlagType);
    int AllocateAndRead(int,int,int);
    bool PositionFileAtFlag(AmberParmFlagType);
    bool PositionFileAtFlag(const char*);

    static void CheckNameWidth(const char*, NameType const&);
    static int AmberIfbox(const Box&);
    int WriteFlagAndFormat(const char*, size_t);
    int WriteSetup(AmberParmFlagType,size_t);
    int WriteInteger(AmberParmFlagType,std::vector<int>const&);
    int WriteDouble(AmberParmFlagType,std::vector<double>const&);
    int WriteDoubleArray(std::vector<double>const&);
    int WriteName(AmberParmFlagType,std::vector<NameType>const&);

    size_t GetFortranBufferSize(int,int,int);
    bool SetFortranType();
};
#endif
