#include "Traj_SQM.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString

void Traj_SQM::WriteHelp() {
  mprintf("\tcharge <c>: Set total integer charge. If not specified it will be calculated from"
          " atomic charges.\n");
}

// Traj_SQM::processWriteArgs()
int Traj_SQM::processWriteArgs(ArgList& argIn) {
  if (argIn.Contains("charge")) {
    charge_ = argIn.getKeyInt( "charge", 0 );
    chargeIsSet_ = true;
  } else
    chargeIsSet_ = false;
  return 0;
}

// Traj_SQM::setupTrajout()
int Traj_SQM::setupTrajout(std::string const& fname, Topology* trajParm,
                           CoordinateInfo const& cInfoIn,
                           int NframesToWrite, bool append)
{
  if (trajParm==0) return 1;
  if (append) {
    mprinterr("Error: Append not supported for SQM.\n");
    return 1;
  }
  SetCoordInfo( cInfoIn );
  if (outfile_.SetupWrite( fname, debug_ )) return 1;
  sqmParm_ = trajParm;
  if (NframesToWrite==1) singleWrite_ = true;
  // Set up title
  std::string outTitle = Title();
  if (outTitle.empty()) {
    outTitle.assign("Cpptraj generated SQM input");
  } else {
    if ( outTitle.size() > 80) {
      mprintf("Warning: Amber SQM title for '%s' too long: truncating.\n[%s]\n",
              outfile_.Filename().base(), outTitle.c_str());
      outTitle.resize(80);
    }
  }
  SetTitle( outTitle );
  // If charge not set, try to determine charge.
  // TODO: Warn if not integer charge
  if (!chargeIsSet_) {
    mprintf("Warning: No charge specified; attempting to calculate charge.\n");
    double qtotal = 0.0;
    for (int i = 0; i < sqmParm_->Natom(); i++)
      qtotal += (*sqmParm_)[i].Charge();
    charge_ = (int)qtotal;
  }
  // Set up header
  header_.assign(" &qmmm\n"
                 "  qm_theory='AM1', qmcharge = "+integerToString(charge_)+", maxcyc = 0,\n"
                 "  tight_p_conv = 1, scfconv = 1.0e-10, pseudo_diag = 0, errconv = 1.0e-10\n"
                 " /\n");
  return 0;
}

// Traj_SQM::writeFrame()
int Traj_SQM::writeFrame(int set, Frame const& frameOut) {
  // If just writing 1 frame dont modify output filename
  if (singleWrite_) {
    if (outfile_.OpenFile()) return 1;
  } else {
    if (outfile_.OpenWriteNumbered( set + 1 )) return 1;
  }
  outfile_.Printf("%s\n", Title().c_str());
  outfile_.Printf("%s", header_.c_str());
  const Topology& parm = static_cast<const Topology&>( *sqmParm_ );
  for (int atom = 0; atom < parm.Natom(); atom++) {
    const double* XYZ = frameOut.XYZ( atom );
    outfile_.Printf("%2d %-4s %12.7f %12.7f %12.7f\n", parm[atom].AtomicNumber(),
                    parm[atom].c_str(), XYZ[0], XYZ[1], XYZ[2]);
  }
  outfile_.CloseFile();
  return 0;
}

// Traj_SQM::Info()
void Traj_SQM::Info() {
  mprintf("is an SQM input file");
}
