// ProgressBar
#include "ProgressBar.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
ProgressBar::ProgressBar() :
  unknown_(0),
  max_(0),
  C_over_max_(1.0),
  targetPercent_(0.0),
  unknownframes_(true)
{}

void ProgressBar::SetupProgress(int maxIn) {
  max_ = maxIn - 1;
  C_over_max_ = 1.0;
  targetPercent_ = 0.0;
  unknownframes_ = false;
  if (max_ < 0) {
    unknownframes_ = true;
    if (maxIn == -1)
      unknown_ = 200; // default
    else
      unknown_ = -maxIn;
    mprintf("\tProgress: '+' = %i iterations.\n", unknown_);
    targetPercent_ = -1.0;
    max_ = unknown_ * 25;
  } else {
    if ( max_ == 0) // One frame
      C_over_max_ = 100.0;
    else 
      C_over_max_ = 100.0 / (float)max_;
  }
}

// CONSTRUCTOR
/** Takes value needed to be considered complete as input. If expected
  * max is less than 0 assume max number of iterations is unknown.
  */
ProgressBar::ProgressBar(int maxIn) {
  SetupProgress( maxIn );
}

// ProgressBar::Update()
/** If current percent is greater than target percent, print current
  * percent complete. If max # is unknown just print '+'.
  */
void ProgressBar::Update(int current) {
  if (unknownframes_) {
    if (targetPercent_ < 0) {
      mprintf("%10i ",current);
      mflush();
      targetPercent_ = (float)(unknown_ - 1);
    }
    if ((float)current > targetPercent_) {
      mprintf("+");
      int target = (int)targetPercent_ + 1;
      if ((target % max_) == 0)
        mprintf("\n%10i ",current);
      targetPercent_ += (float)unknown_;
      mflush();
    } 
  } else {
    if ( current >= max_ )
      mprintf("100%% Complete.\n");
    else {
      float currentPercent = (float)current * C_over_max_;
      if (currentPercent >= targetPercent_) {
        mprintf("%2.0f%% ", currentPercent);
        targetPercent_ += 10.0;
        mflush();
      }
    }
  }
}
// -----------------------------------------------------------------------------
ParallelProgress::ParallelProgress(const ParallelProgress& rhs) :
  C_over_max_(rhs.C_over_max_), tgt_(rhs.tgt_), thread_(rhs.thread_) {}

void ParallelProgress::printProgress(int it) {
  float currentPercent = (float)it * C_over_max_;
  if (currentPercent >= tgt_) {
    mprintf("%2.0f%% ", currentPercent);
    mflush();
    //while (tgt_ <= currentPercent) 
      tgt_ += 10.0;
  }
}

void ParallelProgress::Finish() {
  mprintf("Complete.\n");
}
