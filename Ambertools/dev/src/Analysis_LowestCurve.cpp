#include "Analysis_LowestCurve.h"
#include "CpptrajStdio.h"

Analysis_LowestCurve::Analysis_LowestCurve() : points_(0), step_(0.0) {}

void Analysis_LowestCurve::Help() {
  mprintf("\tpoints <# lowest> [step <stepsize>] <dset0> [<dset1> ...]\n"
          "\t[out <file>] [name <setname>]\n"
          "  Calculate a curve of the average of the # lowest points in bins of stepsize.\n");
}

// Analysis_LowestCurve::Setup()
Analysis::RetType Analysis_LowestCurve::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  DataFile* outfile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  points_ = analyzeArgs.getKeyInt("points", -1);
  if (points_ < 1) {
    mprinterr("Error: 'points' must be specified and > 0\n");
    return Analysis::ERR;
  }
  step_ = analyzeArgs.getKeyDouble("step", 1.0);
  std::string setname = analyzeArgs.GetStringKey("name");
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }
  // Create output data sets
  for (Array1D::const_iterator DS = input_dsets_.begin(); DS != input_dsets_.end(); ++DS) {
    DataSet* dsout = datasetlist->AddSetIdx(DataSet::DOUBLE, setname, DS - input_dsets_.begin());
    if (dsout == 0)
      return Analysis::ERR;
    dsout->SetLegend("LC(" + (*DS)->Legend() + ")" );
    output_sets_.push_back( dsout );
    if (outfile != 0) outfile->AddSet( dsout );
  }

  mprintf("    LOWESTCURVE: Calculating curve of average of %i lowest points in bins of size %g.\n",
          points_, step_);
  mprintf("\t%i data sets.\n", input_dsets_.size());
  if (outfile != 0)
    mprintf("\tWriting results to %s\n", outfile->DataFilename().full());
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->Legend().c_str());

  return Analysis::OK;
}

// Analysis_LowestCurve::Analyze()
Analysis::RetType Analysis_LowestCurve::Analyze() {
  typedef std::vector< std::set<double> > Larray;
  OutArray::const_iterator OUT = output_sets_.begin();
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end();
                             ++DS, ++OUT)
  {
    // Determine an appropriate dimension based on the step size.
    Dimension Xdim;
    Xdim.SetMin( (*DS)->Xcrd(0) );
    Xdim.SetMax( (*DS)->Xcrd((*DS)->Size()-1) );
    Xdim.SetStep( step_ );
    if (Xdim.CalcBinsOrStep()) continue;
    Larray bins_( Xdim.Bins() + 1 );
    mprintf("\tSet '%s' has %i bins (%g < %g, %g)\n", (*DS)->Legend().c_str(),
            Xdim.Bins(), Xdim.Min(), Xdim.Max(), Xdim.Step());
    // Bin each data point. Only save points_ lowest points.
    for (unsigned int n = 0; n != (*DS)->Size(); ++n)
    {
      double xval = (*DS)->Xcrd( n );
      long int idx = (long int)((xval - Xdim.Min()) / Xdim.Step());
      //mprintf("DEBUG: xval=%g, bin=%li\n", xval, idx);
      if ( (int)bins_[idx].size() < points_ )
        // Not enough points yet. Just add this one.
        bins_[idx].insert( (*DS)->Dval( n ) );
      else {
        double dval = (*DS)->Dval( n );
        // Already have enough points. Only add if less than highest point.
        std::set<double>::iterator highest = bins_[idx].end(); // FIXME: Use rbegin?
        --highest; // Position at final element
        if (dval < *highest) {
          bins_[idx].erase( highest );
          bins_[idx].insert( dval );
        }
      }
    }
    // Each bin should now contain only the lowest points_ frames.
    // Plot the average of each bin
    for (int idx = 0; idx != Xdim.Bins(); idx++)
    {
      double avg = 0.0;
      if (!bins_[idx].empty()) {
        for (std::set<double>::const_iterator it = bins_[idx].begin();
                                              it != bins_[idx].end();
                                            ++it)
          avg += *it;
        avg /= (double)bins_[idx].size();
      }
      if ((int)bins_[idx].size() < points_)
        mprintf("Warning: For set '%s'; bin %i had less than %i points.\n",
                (*DS)->Legend().c_str(), idx, points_);
      (*OUT)->Add( idx, &avg );
    }
    (*OUT)->SetDim(Dimension::X, Xdim);
  }
  return Analysis::OK;
}
