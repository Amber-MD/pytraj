#include "Analysis_Average.h"
#include "CpptrajStdio.h"

Analysis_Average::Analysis_Average() {}

void Analysis_Average::Help() {
  mprintf("\t<dset0> [<dset1> ...] [out <file>]\n"
          "  Calculate the average, standard deviation, min, and max of given data sets.\n");
}

// Analysis_Average::Setup()
Analysis::RetType Analysis_Average::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  std::string outname = analyzeArgs.GetStringKey("out");
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  mprintf("    AVERAGE: Calculating average of %i data sets.\n",
          input_dsets_.size());
  if (!outname.empty())
    mprintf("\tWriting results to %s\n", outname.c_str());
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->Legend().c_str());
  if (outfile_.OpenWrite( outname )) return Analysis::ERR;

  return Analysis::OK;
}

// Analysis_Average::Analyze()
Analysis::RetType Analysis_Average::Analyze() {
  double stdev, avg;
  outfile_.Printf("#SetNum\tAverage\tStdev\tMin\tMax\tName\n");
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end(); ++DS)
  {
    if ( (*DS)->Size() < 1)
      mprintf("Warning: Set \"%s\" has no data.\n", (*DS)->Legend().c_str());
    else {
      avg = (*DS)->Avg( stdev );
      double min = (*DS)->Min();
      double max = (*DS)->Max();
      outfile_.Printf("%u\t%f\t%f\t%f\t%f\t\"%s\"\n", DS - input_dsets_.begin(),
                       avg, stdev, min, max, (*DS)->Legend().c_str());
    }
  }
  return Analysis::OK;
}
