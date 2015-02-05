#ifndef INC_TRAJOUT_H
#define INC_TRAJOUT_H
#include "TrajectoryFile.h"
#include "Range.h"
#include "ActionFrameCounter.h"
/// Output trajectory class.
// FIXME: InitTrajWrite should also take # frames to write?
class Trajout : public TrajectoryFile {
  public:
    Trajout();
    ~Trajout();
    inline int InitTrajWrite(std::string const&, ArgList&, Topology*,
                             TrajectoryFile::TrajFormatType);
    inline int InitTrajWrite(std::string const&, Topology*, TrajectoryFile::TrajFormatType);
    int InitStdoutTrajWrite(ArgList&, Topology*, TrajectoryFile::TrajFormatType);
    int InitEnsembleTrajWrite(std::string const&, ArgList const&,
                              Topology*, TrajFormatType, int);
    int InitTrajWriteWithArgs(std::string const&, const char*, Topology*,
                               TrajectoryFile::TrajFormatType);
    void EndTraj();
    int WriteFrame(int, Topology*, Frame const&);
    void PrintInfo(int) const;
    bool TrajIsOpen()        const { return trajIsOpen_;         }
    int NumFramesProcessed() const { return numFramesProcessed_; }
  private:
    int InitTrajWrite(std::string const&, ArgList*, Topology*, TrajectoryFile::TrajFormatType);
    int InitTrajout(std::string const&, ArgList*, Topology*, TrajectoryFile::TrajFormatType);

    int numFramesProcessed_;
    TrajectoryIO* trajio_;
    bool trajIsOpen_;                  ///< If true trajectory has been opened.
    bool nobox_;                       ///< If true do not put box information in output traj
    bool append_;                      ///< If true, append to this file.
    bool hasRange_;                    ///< If true a frame range is defined.
    Range FrameRange_;                 ///< List of frame numbers to write.
    Range::const_iterator rangeframe_; ///< If frame range defined, this is next frame in range.
    ActionFrameCounter frameCount_;    ///< Hold start/stop/offset values
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
int Trajout::InitTrajWrite(std::string const& n, ArgList& a, Topology* p,
                           TrajectoryFile::TrajFormatType t)
{
  return InitTrajWrite( n, &a, p, t );
}
int Trajout::InitTrajWrite(std::string const& n, Topology* p, TrajectoryFile::TrajFormatType t)
{
  return InitTrajWrite( n, 0, p, t );
}
#endif
