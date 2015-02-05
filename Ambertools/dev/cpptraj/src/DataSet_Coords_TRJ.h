#ifndef INC_DATASET_COORDS_TRJ_H
#define INC_DATASET_COORDS_TRJ_H
#include "DataSet_Coords.h"
#include "Trajin.h"
/// Used to read frames from disk.
/** This class will essentially be a copy of TrajinList except it can be
  * used to access frames randomly instead of sequentially.
  */
class DataSet_Coords_TRJ : public DataSet_Coords {
  public:
    DataSet_Coords_TRJ();
    ~DataSet_Coords_TRJ();
    static DataSet* Alloc() { return (DataSet*)new DataSet_Coords_TRJ(); }
    int AddSingleTrajin(std::string const&, ArgList&, Topology*);
    int AddInputTraj(Trajin*);
    // ---- DataSet functions -------------------
    size_t Size() const { return maxFrames_;     }
    int Sync()          { return 1;              }
    void Info() const;
    void Add( size_t, const void* )              { return;     }
    // ----- DataSet_Coords functions ------------
    int AllocateCoords(size_t)                   { return 0;   }
    /// Add a frame.
    void AddFrame(Frame const& fIn) { }
    /// Set CRD at position with frame.
    void SetCRD(int idx, Frame const& fIn) { }
    // -------------------------------------------
    /// Get a frame at position.
    void GetFrame(int idx, Frame& fIn);
    /// Get a frame at position corresponding to mask.
    void GetFrame(int idx, Frame& fIn, AtomMask const& mIn);
   private:
      int SetTrjTopology( Topology const& );
      int UpdateTrjFrames(int);

      typedef std::vector<Trajin*> ListType;
      ListType trajinList_; ///< Input trajectories
      Trajin* Traj_;        ///< Current input trajectory. 
      int currentTrajNum_;  ///< # of currently open input trajectory
      int globalOffset_;    ///< Internal offset for converting global index to traj index
      int maxFrames_;       ///< Total read frames
      bool deleteTrajectories_;
      Frame readFrame_;     ///< For reading in with mask
};
#endif
