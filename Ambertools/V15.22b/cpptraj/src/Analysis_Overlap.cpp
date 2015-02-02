#include <cmath> // fabs
#include "Constants.h" // SMALL
#include "Analysis_Overlap.h"
#include "CpptrajStdio.h"
#include "DataSet_1D.h"

Analysis_Overlap::Analysis_Overlap() : ds1_(0), ds2_(0), useDeviation_(false) {}

void Analysis_Overlap::Help() {
  mprintf("\tds1 <ds1> ds2 <ds2> [rmsd]\n");
}

static inline bool check_type(DataSet* ds, int n_ds) {
  if (ds == 0) {
    mprinterr("Error: Data set ds%i not found.\n", n_ds);
    return true;
  }
  if (ds->Type() != DataSet::FLOAT &&
      ds->Type() != DataSet::DOUBLE &&
      ds->Type() != DataSet::INTEGER) {
    mprinterr("Error: %s: bad set type for overlap.\n", ds->Legend().c_str());
    return true;
  }
  return false;
}

Analysis::RetType Analysis_Overlap::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Keywords
  ds1_ = datasetlist->GetDataSet( analyzeArgs.GetStringKey("ds1") );
  if (check_type(ds1_,1)) return Analysis::ERR;
  ds2_ = datasetlist->GetDataSet( analyzeArgs.GetStringKey("ds2") );
  if (check_type(ds2_,2)) return Analysis::ERR;
  useDeviation_ = analyzeArgs.hasKey("rmsd");

  mprintf("    OVERLAP: Between %s and %s\n", ds1_->Legend().c_str(),
          ds2_->Legend().c_str());
  if (useDeviation_)
    mprintf("\tCalculating overlap using RMSD.\n");

  return Analysis::OK;
}

Analysis::RetType Analysis_Overlap::Analyze() {
  if (ds1_->Size() < 1 || ds2_->Size() < 1) {
    mprinterr("Error: One or both data sets empty (ds1=%i, ds2=%i)\n",
              ds1_->Size(), ds2_->Size());
    return Analysis::ERR;
  }
  if (ds1_->Size() != ds2_->Size()) {
    mprinterr("Error: Data set sizes do not match (ds1=%i, ds2=%i)\n",
              ds1_->Size(), ds2_->Size());
    return Analysis::ERR;
  }
  DataSet_1D const& D1 = static_cast<DataSet_1D&>( *ds1_ );
  DataSet_1D const& D2 = static_cast<DataSet_1D&>( *ds2_ );
  if (useDeviation_) {
    // Determine max value out of either set
    double max = D1.Dval(0);
    for (unsigned int i = 0; i < D1.Size(); i++) {
      if (D1.Dval(i) > max) max = D1.Dval(i);
      if (D2.Dval(i) > max) max = D2.Dval(i);
    }
    double sum = 0.0;
    for (unsigned int i = 0; i < D1.Size(); i++) {
      double diff = (D1.Dval(i)/max) - (D2.Dval(i)/max);
      sum += (diff * diff);
    }
    sum /= (double)D1.Size();
    sum = sqrt( sum );
    mprintf("\tNormalized RMSD of %s from %s is %f\n", ds1_->Legend().c_str(),
            ds2_->Legend().c_str(), 1.0 - sum);
  } else {  
    int Npoints = 0;
    double sum = 0.0;
    for (unsigned int i = 0; i < D1.Size(); i++) {
      double val1 = D1.Dval(i);
      double val2 = D2.Dval(i);
      if (fabs(val1) < Constants::SMALL && fabs(val2) < Constants::SMALL) {
        // No data in either set, do not process;
        continue;
      }
      double denominator = val1 + val2;
      if (fabs(denominator) < Constants::SMALL) {
        // Complete opposite, no overlap, but process
        ++Npoints;
        continue;
      }
      //mprintf("\t%4i %8.3f %8.3f %8.3f %8.3f\n",Npoints,val1,val2,denominator,(1.0 - ( fabs(val1 - val2) / denominator ))); // DEBUG
      sum += (1.0 - ( fabs(val1 - val2) / denominator ));
      ++Npoints;
    }
    if (Npoints < 1)
      sum = 0.0;
    else
      sum /= (double)Npoints;
    mprintf("\t%i of %i points had no data.\n", ds1_->Size() - Npoints, ds1_->Size());
    mprintf("\tPercent overlap between %s and %s is %f\n", ds1_->Legend().c_str(),
            ds2_->Legend().c_str(), sum);
  }
  return Analysis::OK;
}
