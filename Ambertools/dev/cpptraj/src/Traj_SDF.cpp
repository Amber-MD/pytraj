#include "Traj_SDF.h"
#include "CpptrajStdio.h"

bool Traj_SDF::ID_TrajFormat(CpptrajFile& fileIn) {
  return SDFfile::ID_SDF(fileIn);
}

void Traj_SDF::closeTraj() {
  file_.CloseFile();
}

int Traj_SDF::openTrajin() {
  if (file_.OpenFile()) return 1;
   // Read header
  if (file_.ReadHeader()) return 1;
  return 0;
}

int Traj_SDF::setupTrajin(std::string const& fname, Topology* trajParm) {
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  if (openTrajin()) return TRAJIN_ERR;
  // Check number of atoms
  if (file_.SDF_Natoms() != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in SDF file is %i, but associated\n"
              "Error:  topology '%s' has %i\n", file_.SDF_Natoms(),
              trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  closeTraj();
  // Limit to one frame reads for now
  return 1;
}

int Traj_SDF::readFrame(int set, Frame& frameIn) {
  // Limited to 1 frame read for now
  if (set > 0) {
    mprinterr("Error: SDF currently only supports reading 1 frame.\n");
    return 1;
  }
  // Read atoms
  double* Xptr = frameIn.xAddress();
  for (int i = 0; i < file_.SDF_Natoms(); i++) {
    if ( file_.SDF_XYZ( Xptr ) ) {
      mprinterr("Error: Could not read atoms from SDF file.\n");
      return 1;
    }
    Xptr += 3;
  }
  return 0;
}

int Traj_SDF::setupTrajout(std::string const& fname, Topology* trajParm,
                           CoordinateInfo const& cInfoIn,
                               int NframesToWrite, bool append)
{
  mprinterr("Error: SDF writes not yet implemented.\n");
  return 1;
}

int Traj_SDF::writeFrame(int set, Frame const& frameOut) {
  mprinterr("Error: SDF writes not yet implemented.\n");
  return 1;
}

void Traj_SDF::Info() {
  mprintf("is an SDF file");
}
