#ifndef INC_ANALYSIS_HIST_H
#define INC_ANALYSIS_HIST_H
#include "Analysis.h"
#include "DataSet_1D.h"
#include "TrajectoryFile.h" // traj3d
// Class: Analysis_Hist
/// Create an N-dimensional histogram from N input datasets
class Analysis_Hist : public Analysis {
  public :
    enum NormMode { NO_NORM = 0, NORM_SUM, NORM_INT };
    Analysis_Hist();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Hist(); }
    static void Help();
    Analysis::RetType Setup(DataSet_1D*, std::string const&, int, std::string const&,
                            bool, double, bool, double, double, int, double, NormMode,
                            DataSetList&, DataFileList&);
    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    int CheckDimension(std::string const&, DataSetList*);
    int setupDimension(ArgList&, DataSet_1D const&, size_t&);
    int CalcFreeE();
    int Normalize();
    // ---------------------------------
    long int BinIndicesToIndex(std::vector<int> const&);
    bool IncrementBinIndices(std::vector<int>&, int, bool&);
    void PrintBins();

    DataFile* outfile_;                  ///< Output DataFile.
    DataSet* hist_;                      ///< Histogram data set.
    std::vector<double> Bins_;           ///< Histogram data - double in case free E calculated
    typedef std::vector<long int> OffType;
    OffType binOffsets_;                 ///< Bin offsets for calculating index.
    std::vector<DataSet_1D*> histdata_;  ///< Array of data sets to be binned.
    std::vector<ArgList> dimensionArgs_; ///< Array of args defining histogram dims
    typedef std::vector<Dimension> HdimType;
    HdimType dimensions_;                ///< Histogram dimensions.

    int debug_;                          ///< Debug level
    bool calcFreeE_;                     ///< If true, calc free E from hist populations.
    double Temp_;                        ///< temperature to calc free E at.
    NormMode normalize_;                 ///< Normalize histogram
    bool gnuplot_;                       ///< For internal write only
    bool circular_;                      ///< If true, wrap histogram dimensions.
    bool nativeOut_;                     ///< If true, use built in output routine.
    std::string outfilename_;            ///< Stored in case internal write used (DIM > 3)
    size_t N_dimensions_;                ///< # of histogram dimensions.
    Dimension default_dim_;
    bool minArgSet_;
    bool maxArgSet_;
    bool calcAMD_;
    DataSet_1D* amddata_;
    std::string traj3dName_;
    std::string parmoutName_;
    TrajectoryFile::TrajFormatType traj3dFmt_;
};
#endif
