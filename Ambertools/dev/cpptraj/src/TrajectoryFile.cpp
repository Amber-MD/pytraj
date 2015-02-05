#include "TrajectoryFile.h"
#include "CpptrajStdio.h"
// All TrajectoryIO classes go here
#include "Traj_AmberCoord.h"
#ifdef BINTRAJ
  #include "Traj_AmberNetcdf.h"
  #include "Traj_AmberRestartNC.h"
#endif
#include "Traj_PDBfile.h"
#include "Traj_AmberRestart.h"
#include "Traj_Mol2File.h"
#include "Traj_Conflib.h"
#include "Traj_CharmmDcd.h"
#include "Traj_Binpos.h"
#include "Traj_GmxTrX.h"
#include "Traj_SQM.h"
#include "Traj_CIF.h"
#include "Traj_SDF.h"
#include "Traj_Tinker.h"

// ----- STATIC VARS / ROUTINES ------------------------------------------------ 
// NOTE: Must be in same order as TrajFormatType
const FileTypes::AllocToken TrajectoryFile::TF_AllocArray[] = {
# ifdef BINTRAJ
  { "Amber NetCDF",       Traj_AmberNetcdf::ReadHelp, Traj_AmberNetcdf::WriteHelp, Traj_AmberNetcdf::Alloc    },
  { "Amber NC Restart",   Traj_AmberRestartNC::ReadHelp, Traj_AmberRestartNC::WriteHelp, Traj_AmberRestartNC::Alloc },
# else
  { "Amber NetCDF",       0, 0, 0                          },
  { "Amber NC Restart",   0, 0, 0                          },
# endif
  { "PDB",                0, Traj_PDBfile::WriteHelp, Traj_PDBfile::Alloc        },
  { "Mol2",               0, Traj_Mol2File::WriteHelp, Traj_Mol2File::Alloc       },
  { "CIF",                0, 0, Traj_CIF::Alloc            },
  { "Charmm DCD",         0, Traj_CharmmDcd::WriteHelp, Traj_CharmmDcd::Alloc      },
  { "Gromacs TRX",        0, 0, Traj_GmxTrX::Alloc         },
  { "BINPOS",             0, 0, Traj_Binpos::Alloc         },
  { "Amber Restart",      Traj_AmberRestart::ReadHelp, Traj_AmberRestart::WriteHelp, Traj_AmberRestart::Alloc   },
  { "Tinker file",        0, 0, Traj_Tinker::Alloc         },
  { "Amber Trajectory",   0, Traj_AmberCoord::WriteHelp, Traj_AmberCoord::Alloc     },
  { "SQM Input",          0, Traj_SQM::WriteHelp, Traj_SQM::Alloc            },
  { "SDF",                0, 0, Traj_SDF::Alloc            },
  { "LMOD conflib",       0, 0, Traj_Conflib::Alloc        },
  { "Unknown trajectory", 0, 0, 0                          }
};

const FileTypes::KeyToken TrajectoryFile::TF_KeyArray[] = {
  { AMBERNETCDF,    "netcdf",    ".nc"      },
  { AMBERNETCDF,    "cdf",       ".nc"      },
  { AMBERRESTARTNC, "ncrestart", ".ncrst"   },
  { AMBERRESTARTNC, "restartnc", ".ncrst"   },
  { PDBFILE,        "pdb",       ".pdb"     },
  { MOL2FILE,       "mol2",      ".mol2"    },
  { CIF,            "cif",       ".cif"     },
  { CHARMMDCD,      "dcd",       ".dcd"     },
  { CHARMMDCD,      "charmm",    ".dcd"     },
  { GMXTRX,         "trr",       ".trr"     },
  { BINPOS,         "binpos",    ".binpos"  },
  { AMBERRESTART,   "restart",   ".rst7"    },
  { AMBERRESTART,   "restrt",    ".rst7"    },
  { AMBERRESTART,   "rest",      ".rst7"    },
  { AMBERTRAJ,      "crd",       ".crd"     },
  { CONFLIB,        "conflib",   ".conflib" },
  { SQM,            "sqm",       ".sqm"     },
  { SDF,            "sdf",       ".sdf"     },
  { UNKNOWN_TRAJ,   0,           0          }
};
// -----------------------------------------------------------------------------

// CONSTRUCTOR
TrajectoryFile::TrajectoryFile() : debug_(0), trajParm_(0) {}

// TrajectoryFile::SetDebug()
/** Set debug level. */
void TrajectoryFile::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("\tTrajectoryFile debug level set to %i\n",debug_);
}

// TrajectoryFile::SetTrajFileName()
void TrajectoryFile::SetTrajFileName(std::string const& full, bool isRead) {
  if (isRead)
    trajName_.SetFileNameWithExpansion( full );
  else
    trajName_.SetFileName( full );
} 

int TrajectoryFile::SetTrajParm( Topology* tparmIn ) {
  if (tparmIn==0) {
    mprinterr("Error: TrajectoryFile: Parm file is null.\n");
    return 1;
  }
  trajParm_ = tparmIn;
  return 0;
}

// TrajectoryFile::DetectFormat()
TrajectoryIO* TrajectoryFile::DetectFormat(std::string const& fname, TrajFormatType& ttype) {
  CpptrajFile file;
  if (file.SetupRead(fname, 0) == 0) {
    for (int i = 0; i < (int)UNKNOWN_TRAJ; i++) {
      ttype = (TrajFormatType)i;
      TrajectoryIO* IO = (TrajectoryIO*)FileTypes::AllocIO(TF_AllocArray, ttype, true );
      if (IO != 0 && IO->ID_TrajFormat( file ))
        return IO;
      delete IO;
    }
  }
  ttype = UNKNOWN_TRAJ;
  return 0;
}
