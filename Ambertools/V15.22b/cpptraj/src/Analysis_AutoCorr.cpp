#include "Analysis_AutoCorr.h"
#include "CpptrajStdio.h"
#include "DataSet_Vector.h"

// CONSTRUCTOR
Analysis_AutoCorr::Analysis_AutoCorr() :
  lagmax_(-1),
  usefft_(true),
  calc_covar_(true)
{}

void Analysis_AutoCorr::Help() {
  mprintf("\t[name <dsetname>] <dsetarg0> [<dsetarg1> ...] [out <filename>]\n"
          "\t[lagmax <lag>] [nocovar] [direct]\n"
          "  Calculate autocorrelation functions for selected data set(s)\n");
}

Analysis::RetType Analysis_AutoCorr::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  const char* calctype;

  std::string setname = analyzeArgs.GetStringKey("name");
  DataFile* outfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  lagmax_ = analyzeArgs.getKeyInt("lagmax",-1);
  calc_covar_ = !analyzeArgs.hasKey("nocovar");
  usefft_ = !analyzeArgs.hasKey("direct");
  // Select datasets from remaining args
  ArgList dsetArgs = analyzeArgs.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa)
    dsets_ += datasetlist->GetMultipleSets( *dsa );
  if (dsets_.empty()) {
    mprinterr("Error: autocorr: No data sets selected.\n");
    return Analysis::ERR;
  }
  // If setname is empty generate a default name
  if (setname.empty())
    setname = datasetlist->GenerateDefaultName( "autocorr" );
  // Setup output datasets
  int idx = 0;
  for (DataSetList::const_iterator DS = dsets_.begin(); DS != dsets_.end(); ++DS) {
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::DOUBLE, setname, idx++ );
    if (dsout==0) return Analysis::ERR;
    dsout->SetLegend( (*DS)->Legend() );
    outputData_.push_back( dsout );
    // Add set to output file
    if (outfile != 0) outfile->AddSet( outputData_.back() );
  }
 
  if (calc_covar_)
    calctype = "covariance";
  else
    calctype = "correlation";
 
  mprintf("    AUTOCORR: Calculating auto-%s for %i data sets:\n", calctype, dsets_.size());
  dsets_.List();
  if (lagmax_!=-1)
    mprintf("\tLag max= %i\n", lagmax_);
  if ( !setname.empty() )
    mprintf("\tSet name: %s\n", setname.c_str() );
  if ( outfile != 0 )
    mprintf("\tOutfile name: %s\n", outfile->DataFilename().base());
  if (usefft_)
    mprintf("\tUsing FFT to calculate %s.\n", calctype);
  else
    mprintf("\tUsing direct method to calculate %s.\n", calctype);

  return Analysis::OK;
}

Analysis::RetType Analysis_AutoCorr::Analyze() {
  for (unsigned int ids = 0; ids < dsets_.size(); ids++) {
    mprintf("\t\tCalculating AutoCorrelation for set %s\n", dsets_[ids]->Legend().c_str());
    DataSet_1D& Ct = static_cast<DataSet_1D&>( *outputData_[ids] );
    if (dsets_[ids]->Type() == DataSet::VECTOR) {
      DataSet_Vector const& set = static_cast<DataSet_Vector const&>( *dsets_[ids] );
      set.CalcVectorCorr( set, Ct, lagmax_ );
    } else {
      DataSet_1D const& set = static_cast<DataSet_1D const&>( *dsets_[ids] );
      set.CrossCorr( set, Ct, lagmax_, calc_covar_, usefft_ );
    }
  }

  return Analysis::OK;
}
