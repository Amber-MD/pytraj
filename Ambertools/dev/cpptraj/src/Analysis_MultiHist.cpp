#include "Analysis_MultiHist.h"
#include "Analysis_Hist.h"
#include "Analysis_KDE.h"
#include "CpptrajStdio.h"
#include "Array1D.h"

Analysis_MultiHist::Analysis_MultiHist() {}

Analysis_MultiHist::~Analysis_MultiHist() {
  for (Harray::const_iterator ana = Histograms_.begin(); ana != Histograms_.end(); ++ana)
    delete *ana;
}

void Analysis_MultiHist::Help() {
  mprintf("\t[out <filename>] [name <dsname>] [norm | normint] [kde]\n"
          "\t[min <min>] [max <max>] [step <step>] [bins <bins>] [free <T>]\n"
          "\t <dsetarg0> [ <dsetarg1> ... ]\n"
          "  Histogram each data set separately in 1D.\n");
}

Analysis::RetType Analysis_MultiHist::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  bool useKdehist = analyzeArgs.hasKey("kde");
  double min = 0.0;
  double max = 0.0;
  Analysis_Hist::NormMode normalize = Analysis_Hist::NO_NORM;
  std::string outfilename = analyzeArgs.GetStringKey("out");
  std::string setname = analyzeArgs.GetStringKey("name");
  bool minArgSet = analyzeArgs.Contains("min");
  if (minArgSet) min = analyzeArgs.getKeyDouble("min", 0.0);
  bool maxArgSet = analyzeArgs.Contains("max");
  if (maxArgSet) max = analyzeArgs.getKeyDouble("max", 0.0);
  double step = analyzeArgs.getKeyDouble("step", -1.0);
  int bins = analyzeArgs.getKeyInt("bins", -1);
  if (analyzeArgs.hasKey("norm"))
    normalize = Analysis_Hist::NORM_SUM;
  else if (analyzeArgs.hasKey("normint"))
    normalize = Analysis_Hist::NORM_INT;
  double Temp = analyzeArgs.getKeyDouble("free",-1.0);
  // Remaining args should be data sets
  Array1D inputDsets;
  if (inputDsets.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  // Create a histogram analysis for each data set.
  Histograms_.clear();
  for (Array1D::const_iterator ds = inputDsets.begin(); ds != inputDsets.end(); ++ds)
  {
    Analysis::RetType err = Analysis::OK;
    Analysis* ana = 0;
    if (useKdehist) {
      Analysis_KDE* k_ana = new Analysis_KDE();
      err = k_ana->Setup( (*ds), setname, ds - inputDsets.begin(), outfilename, 
                          minArgSet, min, maxArgSet, max, step, bins, Temp,
                          *datasetlist, *DFLin );
      ana = (Analysis*)k_ana;
    } else {
      Analysis_Hist* h_ana = new Analysis_Hist();
      err = h_ana->Setup( (*ds), setname, ds - inputDsets.begin(), outfilename,
                          minArgSet, min, maxArgSet, max, step, bins, Temp, 
                          normalize, *datasetlist, *DFLin );
      ana = (Analysis*)h_ana;
    }
    if (err != Analysis::OK) {
      mprinterr("Error: Could not set up histogram for %s\n", (*ds)->Legend().c_str());
      delete ana;
      return Analysis::ERR;
    }
    Histograms_.push_back( ana );
  }
  if (Histograms_.empty()) {
    mprinterr("Error: No histograms defined.\n");
    return Analysis::ERR;
  }
  mprintf("    MULTIHIST: Creating 1D histograms for %zu data sets:\n\t", Histograms_.size());
  for (Array1D::const_iterator ds = inputDsets.begin(); ds != inputDsets.end(); ++ds)
    mprintf(" %s", (*ds)->Legend().c_str());
  mprintf("\n");
  return Analysis::OK;
}

Analysis::RetType Analysis_MultiHist::Analyze() {
  for (Harray::const_iterator ana = Histograms_.begin(); ana != Histograms_.end(); ++ana)
    (*ana)->Analyze();
  return Analysis::OK;
}
