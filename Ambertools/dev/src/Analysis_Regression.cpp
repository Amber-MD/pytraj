#include "Analysis_Regression.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"

Analysis_Regression::Analysis_Regression() {}

void Analysis_Regression::Help() {
  mprintf("\t<dset0> [<dset1> ...] [name <name>] [out <filename>]\n"
          "  Calculate linear regression lines for given data sets.\n");
}

// Analysis_Regression::Setup()
Analysis::RetType Analysis_Regression::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  DataFile* outfile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
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
  // Setup output data sets
  int idx = 0;
  if ( input_dsets_.size() == 1 )
    idx = -1; // Only one input set, no need to refer to it by index
  // If setname is empty generate a default name
  if (setname.empty())
    setname = datasetlist->GenerateDefaultName( "LR" );
  for ( Array1D::const_iterator DS = input_dsets_.begin();
                                DS != input_dsets_.end(); ++DS)
  {
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::XYMESH, setname, idx++ );
    if (dsout==0) return Analysis::ERR;
    dsout->SetLegend( "LR(" + (*DS)->Legend() + ")" );
    output_dsets_.push_back( (DataSet_1D*)dsout );
    if (outfile != 0) outfile->AddSet( dsout );
  }

  mprintf("    REGRESSION: Calculating linear regression of %i data sets.\n",
          input_dsets_.size());
  if (outfile != 0)
    mprintf("\tOutput to %s\n", outfile->DataFilename().full());
  //if (!outname.empty())
  //  mprintf("\tWriting results to %s\n", outname.c_str());
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->Legend().c_str());
  //if (outfile_.OpenWrite( outname )) return Analysis::ERR;

  return Analysis::OK;
}

// Analysis_Regression::Analyze()
Analysis::RetType Analysis_Regression::Analyze() {
  int nerr = 0;
  //outfile_.Printf("#SetNum\tAverage\tStdev\tMin\tMax\tName\n");
  Array1D::const_iterator dsout = output_dsets_.begin();
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end();
                             ++DS, ++dsout)
  {
    if ( (*DS)->Size() < 2)
      mprintf("Warning: Set \"%s\" does not have enough data for regression (%zu points).\n", 
              (*DS)->Legend().c_str(), (*DS)->Size());
    else {
      DataSet_Mesh mesh;
      double slope, intercept, correl;
      // Set XY mesh
      mesh.SetMeshXY( *(*DS) );
      mprintf("  %zu: %s\n", DS - input_dsets_.begin(), (*DS)->Legend().c_str());
      int err = mesh.LinearRegression( slope, intercept, correl, false );
      nerr += err;
      if (err == 0) {
        // Calculate fitted function
        DataSet_Mesh& outMesh = static_cast<DataSet_Mesh&>( *(*dsout) );
        for (unsigned int i = 0; i < mesh.Size(); i++) {
          double x = mesh.X( i );
          outMesh.AddXY( x, slope * x + intercept );
        }
      }
    }
  }
  if (nerr > 0) return Analysis::ERR;
  return Analysis::OK;
}
