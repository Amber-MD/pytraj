// Traj_AmberRestart
#include <cstdio> // sscanf
#include "Traj_AmberRestart.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NoTrailingWhitespace

// CONSTRUCTOR
Traj_AmberRestart::Traj_AmberRestart() :
  restartAtoms_(0),
  natom3_(0),
  numBoxCoords_(0),
  restartTime_(-1.0),
  restartTemp_(-1.0),
  time0_(1.0),
  dt_(1.0),
  singleWrite_( false),
  readAccess_(false),
  useVelAsCoords_(false)
{}

/** Check for an integer (I5) followed by 0-2 scientific floats (E15.7) */
bool Traj_AmberRestart::ID_TrajFormat(CpptrajFile& fileIn) {
  // Assume file set up for read
  if (fileIn.OpenFile()) return false;
  if (fileIn.NextLine()==0) return false; // Title
  std::string buffer2 = fileIn.GetLine(); // natom, [time, temp]
  fileIn.CloseFile();
  if (buffer2.size() <= 36) {
    int i=0;
    for (; i<5; i++) {
      if (!isspace(buffer2[i]) && !isdigit(buffer2[i])) break;
    }
    if ( i==5 ) {
      if (debug_>0) mprintf("  AMBER RESTART file\n");
      return true;
    }
  }
  return false;
}

// Traj_AmberRestart::closeTraj()
void Traj_AmberRestart::closeTraj() {
  file_.CloseFile();
}

// Traj_AmberRestart::openTrajin()
int Traj_AmberRestart::openTrajin() {
  if (file_.OpenFile()) return 1;
  // Read past title and natom/time/Temp lines
  if (file_.NextLine() == 0) return 1;
  if (file_.NextLine() == 0) return 1;
  return 0; 
}

void Traj_AmberRestart::WriteHelp() {
  mprintf("\tnovelocity: Do not write velocities to restart file.\n"
          "\tremdtraj:   Write temperature to restart file (will also write time).\n"
          "\ttime0:      Time for first frame (if not specified time is not written).\n"
          "\tdt:         Time step for subsequent frames, t=(time0+frame)*dt; (default 1.0)\n");
}

// Traj_AmberRestart::processWriteArgs()
int Traj_AmberRestart::processWriteArgs(ArgList& argIn) {
  // For write, assume we want velocities unless specified
  SetVelocity( !argIn.hasKey("novelocity") );
  SetTemperature(argIn.hasKey("remdtraj"));
  time0_ = argIn.getKeyDouble("time0", -1.0);
  dt_ = argIn.getKeyDouble("dt",1.0);
  // If temperature requested write time as well or format will break.
  if (HasT() && time0_ < 0.0)
    time0_ = 1.0;
  return 0;
}

// Traj_AmberRestart::setupTrajout()
/** Allocate a character buffer based on number of coords and whether 
  * velocities/box info is present.
  */
int Traj_AmberRestart::setupTrajout(std::string const& fname, Topology* trajParm, 
                                    int NframesToWrite, bool append)
{
  if (append) {
    mprinterr("Error: Append not supported for Amber Restart.\n");
    return 1;
  }
  if (file_.SetupWrite( fname, debug_ )) return 1;
  readAccess_ = false;
  // Set trajectory info
  restartAtoms_ = trajParm->Natom();
  natom3_ = restartAtoms_ * 3;
  // Calculate the length of coordinate frame in bytes
  file_.SetupFrameBuffer( natom3_, 12, 6 ); 
  // Dont know ahead of time if velocities will be used, allocate space
  // just in case. Velocity will not be written if V input is null.
  file_.ResizeBuffer( natom3_ );
  // If box coords are present, allocate extra space for them
  if (HasBox()) {
    numBoxCoords_ = 6;
    file_.ResizeBuffer( numBoxCoords_ );
  }
  // If number of frames to write == 1 set singleWrite so we dont append
  // frame # to filename.
  if (NframesToWrite==1) singleWrite_ = true;
  // Set up title
  std::string outTitle = Title();
  if (outTitle.empty()) {
    outTitle.assign("Cpptraj Generated Restart");
    outTitle.resize(80, ' ');
  } else {
    if ( outTitle.size() > 80) {
      mprintf("Warning: Amber restart title for %s too long: truncating.\n[%s]\n",
              file_.Filename().base(), outTitle.c_str());
      outTitle.resize(80);
    }
  }
  SetTitle( outTitle );

  return 0;
}

// Traj_AmberRestart::getBoxAngles()
/** Based on input buffer, determine num box coords and get box angles.
  */
int Traj_AmberRestart::getBoxAngles(std::string const& boxline, Box& trajBox) {
  double box[6];
  if (boxline.empty()) {
    mprinterr("Internal Error: Restart box line is empty.\n");
    return 1;
  }
  numBoxCoords_ = sscanf(boxline.c_str(), "%12lf%12lf%12lf%12lf%12lf%12lf",
                         box, box+1, box+2, box+3, box+4, box+5);
  if (debug_>0) {
    mprintf("DEBUG: Restart BoxLine [%s]\n",boxline.c_str());
    mprintf("       Restart numBoxCoords_=%i\n",numBoxCoords_);
  }
  if (numBoxCoords_==-1) {
    // This can occur if there is an extra newline or whitespace at the end
    // of the restart. Warn the user.
    mprintf("Warning: Restart [%s] appears to have an extra newline or whitespace.\n",
            file_.Filename().base());
    mprintf("         Assuming no box information present.\n");
    trajBox.SetNoBox();
    numBoxCoords_ = 0;
  } else if (numBoxCoords_==6) {
    trajBox.SetBox(box);
  } else {
    mprinterr("Error: Expected 6 box coords in restart box coord line, got %i.\n",
              numBoxCoords_);
    return 1;
  }
  return 0;
}

void Traj_AmberRestart::ReadHelp() {
  mprintf("\tusevelascoords: Use velocities in place of coordinates.\n");
}

int Traj_AmberRestart::processReadArgs(ArgList& argIn) {
  useVelAsCoords_ = argIn.hasKey("usevelascoords");
  return 0;
}

// Traj_AmberRestart::setupTrajin()
/** Set up amber restart file for reading. Check that number of atoms matches
  * number of atoms in associated parmtop. Check for box/velocity info.
  */
int Traj_AmberRestart::setupTrajin(std::string const& fname, Topology* trajParm)
{
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  if (file_.OpenFile()) return TRAJIN_ERR;
  readAccess_ = true;
  // Read in title
  std::string title = file_.GetLine();
  SetTitle( NoTrailingWhitespace(title) );
  // Read in natoms, time, and Replica Temp if present
  std::string nextLine = file_.GetLine();
  if (nextLine.empty()) {
    mprinterr("Error: AmberRestart::open(): Reading restart atoms/time.\n");
    return TRAJIN_ERR;
  }
  SetTemperature(false);
  int nread = sscanf(nextLine.c_str(),"%i %lE %lE",&restartAtoms_,&restartTime_,&restartTemp_);
  if (nread < 1) {
    mprinterr("Error: AmberRestart::open(): Getting restart atoms/time.\n");
    return TRAJIN_ERR;
  } else if (nread == 1) {
    restartTime_ = 0.0;
    restartTemp_ = -1.0;
  } else if (nread == 2) {
    restartTemp_ = -1.0;
  } else {
    SetTemperature( true );
  }
  if (debug_ > 0) 
    mprintf("\tAmber restart: Atoms=%i Time=%lf Temp=%lf\n",restartAtoms_,
            restartTime_, restartTemp_);
  // Check that natoms matches parm natoms
  if (restartAtoms_ != trajParm->Natom()) {
    mprinterr("Error: Number of atoms in Amber Restart %s (%i) does not\n",
              file_.Filename().base(), restartAtoms_);
    mprinterr("       match number in associated parmtop (%i)\n",trajParm->Natom());
    return TRAJIN_ERR;
  }
  natom3_ = restartAtoms_ * 3;
  // Calculate the length of coordinate frame in bytes
  file_.SetupFrameBuffer( natom3_, 12, 6 );
  coordSize_ = file_.FrameSize();
  // Read past restart coords 
  if ( file_.ReadFrame() ) {
    mprinterr("Error: AmberRestart::setupTrajin(): Error reading coordinates.\n");
    return TRAJIN_ERR; 
  }
  // Attempt a second read to get velocities or box coords
  Box boxInfo;
  nread = file_.AttemptReadFrame();
  if ( nread < 0 ) {
    mprinterr("Error: Error attempting to read box line of Amber restart file.\n");
    return TRAJIN_ERR;
  }
  size_t readSize = (size_t)nread;
  //mprintf("DEBUG: Restart readSize on second read = %i\n",readSize);
  if (readSize == 0) {
    // If 0 no box or velo 
    SetVelocity( false );
  } else if (readSize == file_.FrameSize()) {
    // If filled framebuffer again, has velocity info. 
    SetVelocity(true);
    // If we can read 1 more line after velocity, should be box info.
    nextLine = file_.GetLine();
    if (!nextLine.empty()) {
      if (getBoxAngles(nextLine, boxInfo)) return TRAJIN_ERR;
    } 
  } else if (readSize<82) {
    // If we read something but didnt fill framebuffer, should have box coords.
    SetVelocity( false );
    nextLine.assign(file_.Buffer(), readSize);
    if (getBoxAngles(nextLine, boxInfo)) return TRAJIN_ERR;
  } else {
    // Otherwise, who knows what was read?
    mprinterr("Error: AmberRestart::setupTrajin(): When attempting to read in\n");
    mprinterr("Error: box coords/velocity info got %lu chars, expected 0, 37,\n",readSize);
    mprinterr("Error: 73, or %lu.\n",file_.FrameSize());
    mprinterr("Error: This usually indicates a malformed or corrupted restart file.\n");
    return TRAJIN_ERR;
  }
  SetBox( boxInfo );
  if (useVelAsCoords_ && !HasV()) {
    mprinterr("Error: 'usevelascoords' specified but no velocities in this restart.\n");
    return TRAJIN_ERR;
  }
  // Recalculate the frame size
  if (HasV())
    file_.ResizeBuffer( natom3_ );
  if (boxInfo.HasBox())
    file_.ResizeBuffer( numBoxCoords_ );
  file_.CloseFile();
  // Only 1 frame in restart by definition
  return 1;
}

// Traj_AmberRestart::readFrame()
/** Get the restart file frame. If velocities are present, read those too.
  */
int Traj_AmberRestart::readFrame(int set, Frame& frameIn) {
  // Read restart coords into frameBuffer_
  if ( file_.ReadFrame() ) {
    mprinterr("Error: AmberRestart::readFrame(): Error reading coordinates.\n");
    return 1;
  }
  // Set frame temp
  if (HasT())
    *(frameIn.tAddress()) = restartTemp_;
  // Get coords from buffer
  file_.BufferBegin();
  file_.BufferToDouble(frameIn.xAddress(), natom3_);
  // Get velocity from buffer if present
  if (HasV()) {
    if (frameIn.HasVelocity()) {
      if (useVelAsCoords_) 
        file_.BufferToDouble(frameIn.xAddress(), natom3_);
      else
        file_.BufferToDouble(frameIn.vAddress(), natom3_);
    } else
      file_.AdvanceBuffer( coordSize_ );
  }
  // Get box from buffer if present
  if (numBoxCoords_!=0) 
    file_.BufferToDouble(frameIn.bAddress(), numBoxCoords_);

  return 0;
}

// Traj_AmberRestart::readVelocity()
int Traj_AmberRestart::readVelocity(int set, Frame& frameIn) {
  if (HasV()) {
    if ( file_.ReadFrame() ) {
      mprinterr("Error: AmberRestart::readVelocity(): Error reading file.\n");
      return 1;
    }
    // Start buffer right after coords.
    file_.BufferBeginAt(coordSize_);
    file_.BufferToDouble(frameIn.vAddress(), natom3_);
    return 0;
  }
  return 1;
}

// Traj_AmberRestart::writeFrame()
/** Write coords in Frame to file in amber restart format. */
int Traj_AmberRestart::writeFrame(int set, Frame const& frameOut) {
  // If just writing 1 frame dont modify output filename
  if (singleWrite_) {
    if (file_.OpenFile()) return 1;
  } else {
    if (file_.OpenWriteNumbered( set + 1 ) ) return 1;
  }

  // Write out title
  file_.Printf("%-s\n", Title().c_str());
  // Write out atoms
  file_.Printf("%5i",restartAtoms_);
  // Write out restart time
  if (time0_>=0) {
    restartTime_ = (time0_ + (double)set) * dt_;
    file_.Printf("%15.7lE",restartTime_);
  }
  // Write out temperature
  if (HasT())
    file_.Printf("%15.7lE",frameOut.Temperature());
  file_.Printf("\n");

  // Write coords to buffer
  file_.BufferBegin();
  file_.DoubleToBuffer(frameOut.xAddress(), natom3_, "%12.7f");
  // Write velocity to buffer. Check V since velocity not known ahead of time
  if (HasV() && frameOut.HasVelocity())
    file_.DoubleToBuffer(frameOut.vAddress(), natom3_, "%12.7f");
  // Write box to buffer
  if (numBoxCoords_!=0)
    file_.DoubleToBuffer(frameOut.bAddress(), numBoxCoords_, "%12.7f");

  if (file_.WriteFrame()) return 1;

  file_.CloseFile();

  return 0;
}

// Traj_AmberRestart::Info()
void Traj_AmberRestart::Info() {
  mprintf("is an AMBER restart file");
  if (readAccess_) {
    // If read access we know for sure whether there are velocities.
    if (HasV())
      mprintf(" with velocity info");
    else
      mprintf(", no velocities");
    if (useVelAsCoords_) mprintf(" (using velocities as coords)");
  } else {
    // If write, not sure yet whether velocities will be written since
    // it also depends on if the frame has velocity info, so only state
    // if novelocity was specified.
    if (!HasV()) mprintf(", no velocities");
  }
}
