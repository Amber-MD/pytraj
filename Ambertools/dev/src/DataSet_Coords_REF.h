#ifndef INC_DATASET_COORDS_REF_H
#define INC_DATASET_COORDS_REF_H
#include "DataSet_Coords.h"
#include "ArgList.h"
///< Store a single reference frame in double precision.
class DataSet_Coords_REF : public DataSet_Coords {
  public:
    DataSet_Coords_REF() : DataSet_Coords(REF_FRAME), num_(-1) {}
    static DataSet* Alloc() { return (DataSet*) new DataSet_Coords_REF(); }
    // ----- DataSet functions -------------------
    // NOTE: Technically a 1D data set so return 1 if not empty.
    size_t Size() const { if (!frame_.empty()) return 1; else return 0; }
    int Sync()          { return 1;              }
    void Info() const;
    void Add( size_t, const void* )              { return;     }
    // ----- DataSet_Coords functions ------------
    // Size is only ever 1, no need to allocate.
    int AllocateCoords(size_t)                   { return 0;   }
    /// Add a frame.
    inline void AddFrame(Frame const& fIn) { frame_ = fIn; }
    /// Get a frame at position.
    inline void GetFrame(int idx, Frame& fIn) { fIn = frame_; }
    /// Get a frame at position corresponding to mask.
    inline void GetFrame(int idx, Frame& fIn, AtomMask const& mIn) {
      fIn.SetFrame(frame_, mIn);
    }
    /// Set CRD at position with frame.
    inline void SetCRD(int idx, Frame const& fIn) { frame_ = fIn; }
    // -------------------------------------------
    int LoadRef(std::string const&, Topology const&, int);
    int SetupRefFrame(std::string const&, std::string const&, Topology const&, ArgList&, int);
    int StripRef(AtomMask const&);
    Frame const& RefFrame()     const { return frame_; }
    FileName const& FrameName() const { return name_ ; }
    int RefIndex()              const { return num_;   }
  private:
    Frame frame_;       ///< Reference coords.
    //Topology* parm_;    ///< Pointer to associated parm in TopologyList. FIXME: Just copy?
    FileName name_;     ///< Ref structure filename.
    //std::string tag_;   ///< Ref structure optional tag.
    int num_;           ///< Internal reference index #.
};
#endif
