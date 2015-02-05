#include "Analysis_AmdBias.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"

Analysis_AmdBias::Analysis_AmdBias() : 
  ds1_(0), 
  Ethresh_(0.0),
  alpha_(0.0),
  bias_(0)
{}

void Analysis_AmdBias::Help() {
  mprintf("\tds <Edata> ethresh <Ethresh> alpha <alpha> out <filename>\n");
}

Analysis::RetType Analysis_AmdBias::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Keywords
  ds1_ = datasetlist->GetDataSet( analyzeArgs.GetStringKey("ds") );
  if (ds1_ == 0) {
    mprinterr("Error: data set not found ('ds <dsname>')\n");
    return Analysis::ERR;
  }
  if (ds1_->Type() != DataSet::FLOAT &&
      ds1_->Type() != DataSet::DOUBLE &&
      ds1_->Type() != DataSet::INTEGER) {
    mprinterr("Error: %s: bad set type for amdbias.\n", ds1_->Legend().c_str());
    return Analysis::ERR;
  }
  Ethresh_ = analyzeArgs.getKeyDouble("ethresh", -1.0);
  if (Ethresh_ <= 0.0) {
    mprinterr("Error: ethresh must be > 0.0 (%f)\n", Ethresh_);
    return Analysis::ERR;
  }
  alpha_ = analyzeArgs.getKeyDouble("alpha", -1.0);
  if (alpha_ <= 0.0) {
    mprinterr("Error: alpha must be > 0.0 (%f)\n", alpha_);
    return Analysis::ERR;
  }
  DataFile* outfile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Set up bias data set
  bias_ = datasetlist->AddSet(DataSet::DOUBLE, analyzeArgs.GetStringNext(), "bias");
  if (bias_ == 0) return Analysis::ERR;
  if (outfile != 0) outfile->AddSet( bias_ );

  mprintf("    AMDBIAS: Using energy in data set %s, ethresh=%.4f, alpha=%.4f\n", 
          ds1_->Legend().c_str(), Ethresh_, alpha_);
  if (outfile != 0)
    mprintf("\tBias energy will be written to %s\n", outfile->DataFilename().base());

  return Analysis::OK;
}

Analysis::RetType Analysis_AmdBias::Analyze() {
  if (ds1_->Size() < 1) {
    mprinterr("Error: Data set is empty\n");
    return Analysis::ERR;
  }
  DataSet_1D const& ds = static_cast<DataSet_1D&>( *ds1_ );
  DataSet_double& Bias = static_cast<DataSet_double&>( *bias_ );
  Bias.Resize( ds.Size() );
  // Calculate bias according to Eboost = (Ethresh - Ex)^2 / (alpha + (Ethresh - Ex))
  for (unsigned int i = 0; i < ds.Size(); i++) {
    double dval = ds.Dval( i );
    if (dval < Ethresh_) {
      double EV = Ethresh_ - dval;
      double Eboost = (EV * EV) / (alpha_ + EV);
      double Ebias = dval + Eboost;
      Bias[i] = Ebias;
    } else
      Bias[i] = dval;
  }
  return Analysis::OK;
}
