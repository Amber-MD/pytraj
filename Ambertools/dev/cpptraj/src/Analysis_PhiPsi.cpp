#include "Analysis_PhiPsi.h"
#include "CpptrajStdio.h"

Analysis_PhiPsi::Analysis_PhiPsi() {}

void Analysis_PhiPsi::Help() {
  mprintf("\t<dsarg0> [<dsarg1> ...] resrange <range> [out <file>]\n"
          "  Calculate the average and standard deviation for phi/psi pairs.\n");
}

// Analysis_PhiPsi::Setup()
Analysis::RetType Analysis_PhiPsi::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  std::string outname = analyzeArgs.GetStringKey("out");
  std::string rangearg = analyzeArgs.GetStringKey("resrange");
  if (rangearg.empty()) {
    mprinterr("Error: Must specify residue number range argument.\n");
    return Analysis::ERR;
  }
  Range range;
  if (range.SetRange(rangearg)) return Analysis::ERR;
  // Treat remaining args as data set args
  DataSetList sets;
  std::string dsarg = analyzeArgs.GetStringNext();
  while (!dsarg.empty()) {
    sets += datasetlist->GetMultipleSets( dsarg );
    dsarg = analyzeArgs.GetStringNext();
  }
  //std::string dsname = analyzeArgs.GetStringKey("name");
  //if (dsname.empty()) {
  //  mprinterr("Error: Must specify data set name for phi/psi data sets.\n");
  //  return Analysis::ERR;
  //}
  if (sets.empty()) {
    mprinterr("Error: No data sets selected. Specify valid data set arguments.\n");
    return Analysis::ERR;
  }
  // Add only phi/psi pairs that are present.
  for (Range::const_iterator res = range.begin(); res != range.end(); ++res)
  {
    DataSet *phi = 0, *psi = 0;
    for (DataSetList::const_iterator ds = sets.begin(); ds != sets.end(); ++ds)
    {
      // TODO: Handle multiple data set names well.
      if ((*ds)->Idx() == *res)
      {
        if ((*ds)->ScalarType() == DataSet::PHI) phi = *ds;
        else if ((*ds)->ScalarType() == DataSet::PSI) psi = *ds;
      }
    }
    if (phi != 0 && psi != 0) {
      input_dsets_.push_back((DataSet_1D*)phi);
      input_dsets_.push_back((DataSet_1D*)psi);
    } else {
      if (phi == 0) mprintf("Warning: PHI not found for residue %i\n", *res);
      if (psi == 0) mprintf("Warning: PSI not found for residue %i\n", *res);
    }
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  mprintf("    PHIPSI: Calculating average/stdev of %i phi/psi pairs (%zu sets).\n",
          input_dsets_.size() / 2, input_dsets_.size());
  if (debugIn > 0)
    for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
      mprintf("\t%s\n", (*set)->Legend().c_str());
  if (!outname.empty())
    mprintf("\tWriting results to %s\n", outname.c_str());
  if (outfile_.OpenWrite( outname )) return Analysis::ERR;

  return Analysis::OK;
}

// Analysis_PhiPsi::Analyze()
Analysis::RetType Analysis_PhiPsi::Analyze() {
  double aphi, apsi, sphi, spsi;
  outfile_.Printf("%-12s %12s %12s %12s %s\n", "#Phi", "Psi", "SD(Phi)", "SD(Psi)", "Legend");
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end(); DS += 2)
  {
    DataSet_1D const& phi = *(*DS);
    DataSet_1D const& psi = *(*(DS+1));
    std::string legend( phi.Legend() + "-" + psi.Legend() );
    //mprintf("\t%s\n", legend.c_str());
    if ( phi.Size() < 1 || psi.Size() < 1)
      mprintf("Warning: Phi/Psi pair \"%s\" has no data.\n", legend.c_str());
    else {
      aphi = phi.Avg( sphi );
      apsi = psi.Avg( spsi );
      outfile_.Printf("%-12.4f %12.4f %12.4f %12.4f \"%s\"\n", 
                       aphi, apsi, sphi, spsi, legend.c_str());
    }
  }
  return Analysis::OK;
}
