#ifndef INC_ANALYSIS_TIMECORR_H
#define INC_ANALYSIS_TIMECORR_H
#include "Analysis.h"
#include "DataSet_Vector.h"
#include "Corr.h"
/** \author Original Code by Alrun N. Koller & H. Gohlke
  * \author Adapted by DRR
  */
class Analysis_Timecorr : public Analysis {
  public:
    Analysis_Timecorr();

    static DispatchObject* Alloc() { return (DispatchObject*)new Analysis_Timecorr(); }
    static void Help();

    Analysis::RetType Setup(ArgList&,DataSetList*,TopologyList*,DataFileList*,int);
    Analysis::RetType Analyze();
  private:
    struct AvgResults {
      double avgr_;
      double rave_;
      double r3iave_;
      double r6iave_;
    };
    enum ModeType { AUTOCORR = 0, CROSSCORR };
    static const char* ModeString[];
    enum DsetOutType { DPLR_R = 0, DPLR_RRIG, DPLR_R3, DPLR_R6, DPLR_NAME,
                       TC_C,       TC_P,      TC_R3R3,   NDSETOUT };
    typedef std::vector<DataSet*> DSarray;
    struct DStoken {
      const char* Aspect;
      const char* Legend;
      DataSet::DataType Type;
    };
    static DStoken Tokens[];

    std::vector<double> CalculateAverages(DataSet_Vector const&, AvgResults&);
    void CalcCorr(int);
    void Normalize( DataSet*, int, double );

    double tstep_;
    double tcorr_;
    int order_;
    ModeType mode_;
    bool dplr_;
    bool norm_;
    bool drct_;
    ComplexArray data1_;
    ComplexArray data2_;
    DataSet_Vector* vinfo1_;
    DataSet_Vector* vinfo2_;
    DSarray DSOut_;
    std::string filename_;
    std::string Plegend_;
    CorrF_FFT pubfft_;
    CorrF_Direct corfdir_;
};
#endif
