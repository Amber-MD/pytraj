#ifndef INC_DATASET_REMLOG_H
#define INC_DATASET_REMLOG_H
#include <map>
#include <vector>
#include "DataSet.h"
/** Store data from Amber REMD log. For each exchange, the replica index and
  * coordinates index for each replica is stored. These indices start from 1.
  */
class DataSet_RemLog : public DataSet {
  public:
    DataSet_RemLog();
    static DataSet* Alloc() { return (DataSet*)new DataSet_RemLog();}
    /// T-REMD temperature map.
    typedef std::map<double,int> TmapType;
    /// Hold info for a single replica at one exchange.
    class ReplicaFrame;
    /// Allocate for given # of replicas
    void AllocateReplicas(int);
    /// Add given replica frame to specified ensemble member.
    void AddRepFrame(int rep, ReplicaFrame const& frm) { ensemble_[rep].push_back(frm); }
    /// \return replica frame at exchange in specified ensemble member.
    ReplicaFrame const& RepFrame(int exch, int rep) const { return ensemble_[rep][exch]; }
    /// \return replica frame at last exchange in specified ensemble member.
    ReplicaFrame const& LastRepFrame(int rep) const { return ensemble_[rep].back(); }
    /// \return number of exchanges
    int NumExchange() const;
    /// \return true if ensemble is valid.
    bool ValidEnsemble() const;
    /// Trim last replica frame.
    void TrimLastExchange(); 
    // ----- DataSet routines --------------------
    size_t Size() const { return ensemble_.size(); }
    int Sync()          { return 1;                }
    void Info()   const { return;                  }
    // TODO: Remove
    void Add( size_t, const void* ) { return;      }
  private:
    /// Hold info for all exchanges of a single replica.
    typedef std::vector<ReplicaFrame> ReplicaArray;
    /// Hold info for all exchanges of all replicas.
    typedef std::vector<ReplicaArray> ReplicaEnsemble;
    ReplicaEnsemble ensemble_;
};
// ----- Public Class Definitions ----------------------------------------------
class DataSet_RemLog::ReplicaFrame {
  public:
    ReplicaFrame() : replicaIdx_(-1), partnerIdx_(-1), coordsIdx_(-1),
                     success_(false), temp0_(0.0),
                     PE_x1_(0.0), PE_x2_(0.0) {}
    ReplicaFrame(int r, int p, int c, bool s, double t0, double pe1, double pe2) :
      replicaIdx_(r), partnerIdx_(p), coordsIdx_(c), success_(s), 
      temp0_(t0), PE_x1_(pe1), PE_x2_(pe2) {}
    int ReplicaIdx() const { return replicaIdx_; }
    int PartnerIdx() const { return partnerIdx_; }
    int CoordsIdx()  const { return coordsIdx_;  }
    bool Success()   const { return success_;    }
    double Temp0()   const { return temp0_;      }
    double PE_X1()   const { return PE_x1_;      }
    double PE_X2()   const { return PE_x2_;      }
  private:
    int replicaIdx_; ///< Position in ensemble.
    int partnerIdx_; ///< Position to exchange to.
    int coordsIdx_;  ///< Coordinate index.
    bool success_;   ///< Successfully exchanged?
    double temp0_;   ///< Replica bath temperature.
    double PE_x1_;   ///< (HREMD) Potential energy with coords 1.
    double PE_x2_;   ///< (HREMD) Potential energy with coords 2.
};
#endif
