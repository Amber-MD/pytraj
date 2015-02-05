#include "Analysis_CrossCorr.h"
#include "DataSet_MatrixFlt.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "DataSet_1D.h"

// CONSTRUCTOR
Analysis_CrossCorr::Analysis_CrossCorr() : outfile_(0), matrix_(0) {}

void Analysis_CrossCorr::Help() {
  mprintf("\t[name <dsetname>] <dsetarg0> [<dsetarg1> ...] [out <filename>]\n"
          "  Calculate matrix of Pearson product-moment correlation\n"
          "  coefficients between selected data sets.\n");
}

// Analysis_CrossCorr::Setup()
Analysis::RetType Analysis_CrossCorr::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  outfile_ = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Select datasets from remaining args
  ArgList dsetArgs = analyzeArgs.RemainingArgs();
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa)
    dsets_ += datasetlist->GetMultipleSets( *dsa );
  if (dsets_.empty()) {
    mprinterr("Error: crosscorr: No data sets selected.\n");
    return Analysis::ERR;
  }
  // Setup output dataset
  matrix_ = datasetlist->AddSet( DataSet::MATRIX_FLT, setname, "crosscorr" );
  if (outfile_ != 0) {
    matrix_->Dim(Dimension::X).SetLabel("DataSets");
    outfile_->AddSet( matrix_ );
  }
  
  mprintf("    CROSSCORR: Calculating correlation between %i data sets:\n", dsets_.size());
  dsets_.List();
  if ( !setname.empty() )
    mprintf("\tSet name: %s\n", setname.c_str() );
  if ( outfile_ != 0 )
    mprintf("\tOutfile name: %s\n", outfile_->DataFilename().base());

  return Analysis::OK;
}

// Analysis_CrossCorr::Analyze()
Analysis::RetType Analysis_CrossCorr::Analyze() {
  DataSet_MatrixFlt& tmatrix = static_cast<DataSet_MatrixFlt&>( *matrix_ );

  int Nsets = dsets_.size();
  mprintf("\tDataSet Legend:\n");
  std::string Ylabels("\"");
  for (int i = 0; i < Nsets; ++i) {
    mprintf("\t\t%8i: %s\n", i+1, dsets_[i]->Legend().c_str());
    //Xlabels_ += (dsets_[i]->Legend() + ",");
    Ylabels += (integerToString(i+1) + ":" + dsets_[i]->Legend() + ",");
  }
  Ylabels += "\"";
  int Nsets1 = Nsets - 1;
  if (tmatrix.AllocateTriangle(Nsets)) return Analysis::ERR;
  for (int i = 0; i < Nsets1; ++i) {
    for (int j = i + 1; j < Nsets; ++j) {
      //mprinterr("DBG:\tCross corr between %i (%s) and %i (%s)\n",
      //          i, dsets_[i]->Legend().c_str(), j, dsets_[j]->Legend().c_str());
      DataSet_1D const& set1 = static_cast<DataSet_1D const&>( *dsets_[i] );
      double corr = set1.CorrCoeff( *((DataSet_1D*)dsets_[j]) );
      tmatrix.AddElement( (float)corr );
    }
  }
  if (outfile_ != 0)
    outfile_->ProcessArgs("ylabels " + Ylabels);

  return Analysis::OK;
}
