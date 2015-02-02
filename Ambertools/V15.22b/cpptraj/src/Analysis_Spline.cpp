#include "Analysis_Spline.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL

Analysis_Spline::Analysis_Spline() :
  outfile_(0),
  meshsize_(0),
  meshmin_(0.0),
  meshmax_(0.0),
  meshfactor_(0.0),
  useDefaultMin_(false),
  useDefaultMax_(false)
{}

void Analysis_Spline::Help() {
  mprintf("\t<dset0> [<dset1> ...] [out <outfile>] [meshsize <n> | meshfactor <x>]\n"
          "\t[meshmin <mmin>] [meshmax <mmax>]\n"
          "  Cubic spline the given data sets.\n");
}

Analysis::RetType Analysis_Spline::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  outfile_ = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  meshsize_ = analyzeArgs.getKeyInt("meshsize", 0);
  meshfactor_ = -1.0;
  if (meshsize_ < 3) {
    meshfactor_ = analyzeArgs.getKeyDouble("meshfactor", -1.0);
    if (meshfactor_ < Constants::SMALL) {
      mprinterr("Error: Either meshsize must be specified and > 2, or meshfactor must be\n"
                "Error:   specified and > 0.0\n");
      return Analysis::ERR;
    }
  }
  if (analyzeArgs.Contains("meshmin")) {
    meshmin_ = analyzeArgs.getKeyDouble("meshmin", 0.0);
    useDefaultMin_ = true;
  } else
    useDefaultMin_ = false;
  if (analyzeArgs.Contains("meshmax")) {
    meshmax_ = analyzeArgs.getKeyDouble("meshmax", -1.0);
    useDefaultMax_ = true;
  } else
    useDefaultMax_ = false;
  if (useDefaultMin_ && useDefaultMax_ && meshmax_ < meshmin_) {
    mprinterr("Error: meshmax must be > meshmin\n");
    return Analysis::ERR;
  }
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  // Set up output datasets
  Dimension Xdim( meshmin_, (meshmax_ - meshmin_) / (double)meshsize_,
                  meshsize_ );
  for (Array1D::const_iterator dsIn = input_dsets_.begin();
                               dsIn != input_dsets_.end(); ++dsIn)
  {
    DataSet* ds = datasetlist->AddSet(DataSet::XYMESH, setname, "Spline");
    if (ds == 0) return Analysis::ERR;
    ds->SetLegend( "Spline(" + (*dsIn)->Legend() + ")" );
    // TODO: Set individually based on input_dsets_
    ds->SetDim(Dimension::X, Xdim);
    if (outfile_ != 0) outfile_->AddSet( ds );
    output_dsets_.push_back( (DataSet_Mesh*)ds );
  }
  /*outfile_->Dim(Dimension::X).SetMin( meshmin_ );
  double meshstep = (meshmax_ - meshmin_) / (double)meshsize_;
  outfile_->Dim(Dimension::X).SetStep( meshstep );*/

  mprintf("    SPLINE: Applying cubic splining to %u data sets\n", input_dsets_.size());
  if (meshfactor_ < 0)
    mprintf("\tMesh size= %i\n", meshsize_);
  else
    mprintf("\tMesh size will be input set size multiplied by %f\n", meshfactor_);
  if (useDefaultMin_)
    mprintf("\tMesh min= %f,", meshmin_);
  else
    mprintf("\tMesh min will be input set min,");
  if (useDefaultMax_)
    mprintf(" Mesh max= %f\n", meshmax_);
  else
    mprintf(" Mesh max will be input set max.\n");
  if (outfile_ != 0) {
    if (!setname.empty())
      mprintf("\tOutput set name: %s\n", setname.c_str());
    mprintf("\tOutfile name: %s\n", outfile_->DataFilename().base());
  }
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->Legend().c_str());
  return Analysis::OK;
}

// TODO: Each data set needs dimension information.
Analysis::RetType Analysis_Spline::Analyze() {
  double mmin, mmax;
  int msize;
  for (unsigned int idx = 0; idx < input_dsets_.size(); idx++) {
    DataSet_1D const& ds = static_cast<DataSet_1D const&>(*input_dsets_[idx]);
    if (useDefaultMin_)
      mmin = meshmin_;
    else
      mmin = ds.Dim(0).Min();
    if (useDefaultMax_)
      mmax = meshmax_;
    else
      mmax = ds.Dim(0).Max();
    if (meshfactor_ > 0)
      msize = (int)((double)ds.Size() * meshfactor_);
    else
      msize = meshsize_;
    // Set up output dimension
    mprintf("\t%s: Setting mesh from %f->%f, size=%i,", ds.Legend().c_str(), mmin, mmax, msize);
    output_dsets_[idx]->SetDim( Dimension::X,
        Dimension(mmin, (mmax - mmin)/(double)msize, msize, ds.Dim(0).Label()) );
    output_dsets_[idx]->CalculateMeshX(msize, mmin, mmax);
    mprintf(" set min=%f, set step=%f\n", ds.Dim(0).Min(), ds.Dim(0).Step());
    if (!ds.Dim(0).MinIsSet())
      mprinterr("Internal Error: %s min is not set!\n", ds.Legend().c_str());
    if (ds.Dim(0).Step() < 0.0)
      mprinterr("Internal Error: %s step is not set!\n", ds.Legend().c_str());
    // Calculate mesh Y values.
    output_dsets_[idx]->SetSplinedMesh( ds );
    // DEBUG
    //for (unsigned int i = 0; i < output_dsets_[idx]->Size(); i++)
    //  mprintf("\t%12.4f %12.4g\n", output_dsets_[idx]->X(i), output_dsets_[idx]->Y(i));
  }
  return Analysis::OK;
}
