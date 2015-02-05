#ifndef INC_DATAIO_REMLOG_H
#define INC_DATAIO_REMLOG_H
#include "DataIO.h"
#include "BufferedLine.h"
#include "DataSet_RemLog.h"
/// Read replica exchange log data.
class DataIO_RemLog : public DataIO {
  public:
    DataIO_RemLog();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_RemLog(); }
    static void ReadHelp();
    int processReadArgs(ArgList&);
    int ReadData(std::string const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&) { return 0; }
    int WriteData(std::string const&, DataSetList const&) { return 1; }
    int WriteData2D(std::string const&, DataSetList const&) { return 1; }
    int WriteData3D(std::string const&, DataSetList const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    // NOTE: Must match ExchgDescription
    enum ExchgType { UNKNOWN = 0, TREMD, HREMD, MREMD, RXSGLD };
    static const char* ExchgDescription[];
    typedef std::vector<std::string> Sarray;
    int ReadRemlogHeader(BufferedLine&, ExchgType&) const;
    int ReadRemdDimFile(std::string const&);
    DataSet_RemLog::TmapType SetupTemperatureMap(BufferedLine&,std::vector<int>&) const;
    int CountHamiltonianReps(BufferedLine&) const;
    int OpenMremdDims(std::vector<BufferedLine>&, Sarray const&);
    void SetupDim1Group( int );
    void PrintReplicaStats(DataSet_RemLog const&);

    Sarray logFilenames_; ///< Replica log file names.
    std::string dimfile_;
    std::string crdidx_;
    int n_mremd_replicas_;
    bool processMREMD_;
    bool searchForLogs_;
    class GroupReplica;
    typedef std::vector<GroupReplica> GroupArray;
    typedef std::vector<GroupArray> GroupDimType;
    std::vector<GroupDimType> GroupDims_;
    std::vector<ExchgType> DimTypes_;
    // Used for getting temps/coord indices from T-remlog
    struct TlogType {
      double t0;
      int crdidx;
    };
    struct TlogType_cmp {
      inline bool operator()(TlogType const& first, TlogType const& second) const {
        return (first.t0 < second.t0);
      }
    };
};
/// Used to hold replica partner info in M-REMD simulations
class DataIO_RemLog::GroupReplica {
  public:
    GroupReplica() : l_partner_(-1), me_(-1), r_partner_(-1) {}
    GroupReplica(const GroupReplica& rhs) :
      l_partner_(rhs.l_partner_), me_(rhs.me_), r_partner_(rhs.r_partner_) {}
    GroupReplica(int l, int m, int r) : l_partner_(l), me_(m), r_partner_(r) {}
    int L_partner() const { return l_partner_; }
    int Me()        const { return me_;        }
    int R_partner() const { return r_partner_; }
  private:
    int l_partner_, me_, r_partner_;
}; 
#endif
