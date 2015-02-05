#include "ActionFrameCounter.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ActionFrameCounter::ActionFrameCounter() : start_(0), stop_(-1), offset_(1) {}

const char* ActionFrameCounter::HelpText = 
  "[start <start>] [stop <stop>] [offset <offset>]";

/** Expected args: [start <start>] [stop <stop>] [offset <offset>] 
  * Defaults: start = 1, stop = -1 (end frame), offset = 1
  */
int ActionFrameCounter::InitFrameCounter(ArgList& argIn) {
  start_ = argIn.getKeyInt("start",1);
  if ( start_ < 1 ) {
    mprintf("Warning: start frame %i is less than 1, setting to 1.\n", start_);
    start_ = 1;
  }
  // User start/stop args begin at frame 1. 
  --start_;
  stop_ = argIn.getKeyInt("stop",-1);
  // For backwards compat., check for 'end' (matrix only?)
  // NOTE: Deprecate?
  if (stop_ == -1)
    stop_ = argIn.getKeyInt("end", -1);
  if (stop_!=-1) {
    --stop_;
    if (stop_ < start_)
      mprintf("Warning: stop frame %i less than start (%i); only 1 frame will be processed.\n",
              stop_ + 1, start_ + 1);
  }
  offset_ = argIn.getKeyInt("offset",1);
  if (offset_ < 1) {
    mprinterr("Error: offset cannot be less than 1 (%i)\n", offset_);
    return 1;
  }
  return 0;
}

void ActionFrameCounter::FrameCounterInfo() const {
  mprintf("\tStart: %i  Stop:", start_ + 1);
  if (stop_ == -1)
    mprintf(" Final frame");
  else
    mprintf(" %i", stop_ + 1);
  if (offset_ > 1)
    mprintf("  Offset: %i", offset_);
  mprintf("\n");
}

void ActionFrameCounter::FrameCounterBrief() const {
  if (stop_ != -1)
    mprintf(" (%i-%i, %i)", start_+1, stop_+1, offset_);
  else
    mprintf(" (%i-Last, %i)", start_+1, offset_);
}
