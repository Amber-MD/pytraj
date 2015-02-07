#include "Trajin.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL

// CONSTRUCTOR
Trajin::Trajin() :
  start_(0),
  stop_(-1),
  offset_(1),
  total_frames_(0),
  total_read_frames_(-1),
  currentFrame_(0),
  numFramesProcessed_(0),
  useProgress_(true),
  isEnsemble_(false)
{}

// Trajin::CheckFrameArgs()
/** Parse argument list for trajectory-related frame args. Frame args start at
  * 1, internal frame #s start at 0. So for a traj with 10 frames:
  * - Internal #: 0 1 2 3 4 5 6 7 8 9
  * - Frame Arg#: 1 2 3 4 5 6 7 8 9 10
  * - Defaults: startArg=1, stopArg=-1, offsetArg=1
  */
int Trajin::CheckFrameArgs(ArgList& argIn, int maxFrames,
                           int& startArg, int& stopArg, int& offsetArg)
{
  if (argIn.hasKey("lastframe")) {
    // lastframe is a special case where only the last frame will be selected
    if (maxFrames > 0) {
      startArg = maxFrames;
      stopArg = maxFrames;
      offsetArg = 1;
    } else {
      mprinterr("Error: lastframe specified but # frames could not be determined.\n");
      return 1;
    }
  } else {
    startArg = argIn.getNextInteger(1);
    // Last explicitly selects final frame as stop arg.
    if (argIn.hasKey("last"))
      stopArg = -1;
    else
      stopArg = argIn.getNextInteger(-1);
    offsetArg = argIn.getNextInteger(1);
  }
  // Check that start argument is valid.
  if (startArg != 1) {
    if (startArg < 1) {
      mprintf("Warning: start argument %i < 1, setting to 1.\n", startArg);
      startArg = 1; //start_ = 0;
    } else if (maxFrames >= 0 && startArg > maxFrames) {
      // startArg==stopArg and greater than # frames, archaic 'lastframe'.
      if (startArg == stopArg) {
        mprintf("Warning: start %i > #Frames (%i), setting to last frame.\n",
                startArg, maxFrames);
        startArg = maxFrames; //start_ = total_frames_ - 1;
      } else {
        mprinterr("Error: start %i > #Frames (%i), no frames will be processed.\n",
                  startArg, maxFrames);
        //start=startArg - 1;
        return 1;
      }
    }
  }
  startArg--; // Internal frame nums start from 0.
  // Check that stop argument is valid
  if (stopArg != -1) {
    if ( (stopArg - 1) < startArg) { // Internal frame nums start from 0.
      mprinterr("Error: stop %i < start, no frames will be processed.\n", stopArg);
      //stop = start;
      return 1;
    } else if (maxFrames >= 0 && stopArg > maxFrames) {
      mprintf("Warning: stop %i > #Frames (%i), setting to max.\n", stopArg, maxFrames);
      stopArg = maxFrames;
    } 
  } else if (maxFrames >= 0) // -1 means use last frame
    stopArg = maxFrames;
  // Check that offset argument is valid.
  if (offsetArg != 1) {
    if (offsetArg < 1) {
      mprintf("Warning: offset %i < 1, setting to 1.\n", offsetArg);
      offsetArg = 1;
    } else if (stopArg != -1 && offsetArg >= (stopArg - startArg)) {
      mprintf("Warning: offset %i is so large that only 1 set will be processed.\n",
              offsetArg);
    }
  }
  //mprintf("DEBUG SetArgs: Start %i Stop %i  Offset %i\n", startArg, stopArg, offsetArg);
  return 0;
}

// Trajin::SetupTrajIO()
int Trajin::SetupTrajIO( std::string const& fname, TrajectoryIO& trajio, ArgList& argIn ) {
  // -1 indicates an error.
  // -2 indicates the number of frames could not be determined, read to EOF.
  total_frames_ = trajio.setupTrajin(fname, TrajParm());
  if (total_frames_ == TrajectoryIO::TRAJIN_ERR) {
    mprinterr("Error: Could not set up %s for reading.\n", fname.c_str());
    return 1;
  }
  if (debug_ > 0) {
    if (total_frames_>-1)
      mprintf("\t'%s' contains %i frames.\n", TrajFilename().base(), total_frames_);
    else
      mprintf("\t'%s' contains an unknown number of frames.\n",TrajFilename().base());
  }
  // Set stop based on calcd number of frames.
  if (total_frames_==0) {
    mprinterr("Error: trajectory %s contains no frames.\n",TrajFilename().base());
    return 1;
  }
  if (total_frames_>0)
    stop_ = total_frames_;
  else
    stop_ = -1;
  // Set the start, stop, and offset args based on user input. Do some bounds
  // checking.
  if (Trajin::CheckFrameArgs( argIn, total_frames_, start_, stop_, offset_))
    return 1;
  return 0;
}

// Trajin::setupFrameInfo()
/** Calculate number of frames that will be read based on start, stop, and
  * offset (total_read_frames). 
  * \return the total number of frames that will be read for this traj.
  * \return -1 if the number of frames could not be determined.
  */
int Trajin::setupFrameInfo() {
  int Nframes;
  int ptraj_start_frame, ptraj_end_frame;
  int traj_start_frame, traj_end_frame;
  // DEBUG - No mpi for now
  int worldrank = 0;
  int worldsize = 1;

  //mprintf("DEBUG: Calling setupFrameInfo for %s with %i %i %i\n",trajName,
  //        start,stop,offset);
  if (stop_==-1) return -1;

  // Calc total frames that will be read
  Nframes = stop_ - start_;
  total_read_frames_ = Nframes / offset_;
  // Round up
  if ( (Nframes % offset_) > 0 )
    ++total_read_frames_;
  //divresult = div( (stop - start), offset);
  //total_read_frames = divresult.quot;
  //if (divresult.rem!=0) total_read_frames++;

  // Calc min num frames read by any given thread
  // last thread gets leftovers
  // In case of 0, last thread gets the frame
  Nframes = total_read_frames_ / worldsize;
  int leftover_frames = total_read_frames_ % worldsize;
  //divresult = div(total_read_frames,worldsize);
  //Nframes=divresult.quot;

  // Ptraj (local) start and end frame
  ptraj_start_frame = (worldrank*Nframes);
  ptraj_end_frame = ptraj_start_frame + Nframes;
  // Last thread gets the leftovers
  if (worldrank == worldsize-1)
    ptraj_end_frame += leftover_frames;

  // Actual Traj start and end frame (for seeking)
  traj_start_frame = (ptraj_start_frame * offset_) + start_;
  traj_end_frame = ((ptraj_end_frame-1) * offset_) + start_;

  start_ = traj_start_frame;
  stop_ = traj_end_frame;

  if ( total_read_frames_ == 0) {
    mprinterr("Error: No frames will be read from %s based on start, stop,\n",TrajFilename().base());
    mprinterr("         and offset values (%i, %i, %i)\n",start_+1,stop_+1,offset_);
  }

  return total_read_frames_;
}

// Trajin::PrepareForRead()
void Trajin::PrepareForRead(bool useIn) {
  numFramesProcessed_ = 0;
  // Setup progress bar
  useProgress_ = useIn;
  if (useProgress_) progress_.SetupProgress( total_read_frames_ );
  currentFrame_ = start_;
}

// Trajin::PrintInfoLine()
void Trajin::PrintInfoLine() const {
  if (stop_ != -1)
    mprintf( "----- %s (%i-%i, %i) -----\n",TrajFilename().base(),start_+1,stop_+1,offset_);
  else
    mprintf( "----- %s (%i-EOF, %i) -----\n",TrajFilename().base(),start_+1,offset_);
}

// Trajin::PrintFrameInfo()
void Trajin::PrintFrameInfo() const {
  if (stop_!=-1 && total_frames_>0)
    //mprintf(": %i-%i, %i (reading %i of %i)",start,stop,offset,total_read_frames,total_frames);
    mprintf(" (reading %i of %i)",total_read_frames_,total_frames_);
  else if (stop_!=-1 && total_frames_ < 0)
    mprintf(" (reading %i)",total_read_frames_);
  else
    mprintf(", unknown #frames, start=%i offset=%i",start_,offset_);
}
