#include "DataIO_VecTraj.h"
#include "CpptrajStdio.h"
#include "DataSet_Vector.h"
#include "ParmFile.h"
#include "Trajout.h"

// CONSTRUCTOR
DataIO_VecTraj::DataIO_VecTraj() : trajoutFmt_(TrajectoryFile::UNKNOWN_TRAJ) {
  SetValid( DataSet::VECTOR );
}

void DataIO_VecTraj::WriteHelp() {
  mprintf("\t[trajfmt <format>] [parmout <file>]\n");
}

int DataIO_VecTraj::processWriteArgs(ArgList& argIn) {
  trajoutFmt_ = TrajectoryFile::GetFormatFromString( argIn.GetStringKey("trajfmt") );
  parmoutName_ = argIn.GetStringKey("parmout");
  return 0;
}

int DataIO_VecTraj::WriteData(std::string const& fname, DataSetList const& SetList) {
  if (SetList.empty()) return 1;
  if (SetList.size() > 1)
    mprintf("Warning: Multiple sets not yet supported for Evecs write.\n");
  DataSet_Vector const& Vec = static_cast<DataSet_Vector const&>( *(*(SetList.begin())) );
  // Create pseudo-topology.
  Topology pseudo;
  pseudo.AddTopAtom(Atom("OXYZ", ' ', 0), 1, "VEC", 0);
  pseudo.AddTopAtom(Atom("VXYZ", ' ', 0), 1, "VEC", 0);
  pseudo.AddBond(0, 1);
  pseudo.CommonSetup(false);
  if (!parmoutName_.empty()) {
    ParmFile pfile;
    if (pfile.WriteTopology( pseudo, parmoutName_, ParmFile::UNKNOWN_PARM, 0 )) {
      mprinterr("Error: Could not write pseudo topology to '%s'\n", parmoutName_.c_str());
      return 1;
    }
  }
  Trajout out;
  if (out.InitTrajWrite(fname, &pseudo, trajoutFmt_) == 0) {
    Frame outFrame(2);
    for (unsigned int i = 0; i != Vec.Size(); ++i) {
      outFrame.ClearAtoms();
      Vec3 const& ovec = Vec.OXYZ(i);
      outFrame.AddVec3( ovec );
      outFrame.AddVec3( Vec[i] + ovec );
      if (out.WriteFrame(i, &pseudo, outFrame)) return 1;
    }
    out.EndTraj();
  } else {
    mprinterr("Error: Could not set up '%s' for write.\n", fname.c_str());
    return 1;
  }
  return 0;
} 
