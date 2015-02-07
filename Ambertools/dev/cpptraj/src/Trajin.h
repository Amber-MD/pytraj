#ifndef INC_TRAJIN_H
#define INC_TRAJIN_H
#include "TrajectoryFile.h"
#include "ProgressBar.h"
/// Class that all input trajectories will inherit.
class Trajin : public TrajectoryFile {
  public:
    Trajin();
    virtual ~Trajin() {}
    virtual int SetupTrajRead(std::string const&, ArgList&, Topology *) = 0;
    virtual int ReadTrajFrame(int, Frame&) = 0;
    virtual int BeginTraj(bool) = 0;
    virtual void EndTraj() = 0;
    virtual void PrintInfo(int) const = 0;
    virtual CoordinateInfo const& TrajCoordInfo() const = 0;

    static int CheckFrameArgs(ArgList&, int, int&, int&, int&);
    inline bool CheckFinished();
    inline void UpdateCounters();
    inline int GetNextFrame(Frame&);
    inline void SetTotalFrames(int); 
    int SetupTrajIO( std::string const&, TrajectoryIO&, ArgList& );
    int setupFrameInfo();
    void PrepareForRead(bool);
    void PrintInfoLine() const;
    void PrintFrameInfo() const;
 
    int TotalFrames()        const { return total_frames_;       }
    int TotalReadFrames()    const { return total_read_frames_;  }
    int CurrentFrame()       const { return currentFrame_;       } // GetNextEnsemble
    int Start()              const { return start_;              }
    int Offset()             const { return offset_;             }
    int NumFramesProcessed() const { return numFramesProcessed_; }
    bool IsEnsemble()        const { return isEnsemble_;         }
    void SetEnsemble(bool b)       { isEnsemble_ = b;            }
  private:
    int start_;              ///< Frame to begin processing
    int stop_;               ///< Frame to end processing
    int offset_;             ///< Number of frames to skip between processed frames
    int total_frames_;       ///< The total number of frames in the traj
    int total_read_frames_;  ///< # frames that will actually be read based on start/stop/offset
    int currentFrame_;       ///< The current frame number being read
    int numFramesProcessed_; ///< Number of frames that have actually been read.
    ProgressBar progress_;   ///< Keep track of trajectory progress
    bool useProgress_;       ///< Indicate whether progress should be shown.
    bool isEnsemble_;        ///< True if this will be processed as an ensemble.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
/** \return true if current frame is out of range. */
bool Trajin::CheckFinished() {
  if (currentFrame_ > stop_ && stop_ != -1) return true;
  if (useProgress_)
    progress_.Update(numFramesProcessed_);
  return false;
}
/** Update current frame and number of frames processed. */
void Trajin::UpdateCounters() {
  ++numFramesProcessed_;
  currentFrame_ += offset_;
}
/** \return 1 if more frames should be read, 0 if finished. */
int Trajin::GetNextFrame(Frame& frameIn) {
  // If the current frame is out of range, exit
  if (CheckFinished()) return 0;
  // Read current frame
  if ( ReadTrajFrame( currentFrame_, frameIn ) ) return 0;
  // Update counters
  UpdateCounters();
  return 1;
}
// Trajin::SetTotalFrames()
void Trajin::SetTotalFrames(int nfIn) { 
  total_frames_ = nfIn;
  if (stop_ > total_frames_) stop_ = total_frames_;
}
#endif
