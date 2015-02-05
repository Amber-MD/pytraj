#include <cmath> // log
#include "Analysis_Divergence.h"
#include "DataSet_1D.h"
#include "CpptrajStdio.h"
#include <limits> // Minimum double val for checking zero

// CONSTRUCTOR
Analysis_Divergence::Analysis_Divergence() : ds1_(0), ds2_(0) {}

void Analysis_Divergence::Help() {
  mprintf("\tds1 <ds1> ds2 <ds2>\n"
          "  Calculate Kullback-Liebler divergence between specified data sets.\n");
}

/// Ensure set is of valid type for divergence calc.
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

// Analysis_Divergence::Setup()
Analysis::RetType Analysis_Divergence::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Keywords
  ds1_ = datasetlist->GetDataSet( analyzeArgs.GetStringKey("ds1") );
  if (check_type(ds1_,1)) return Analysis::ERR;
  ds2_ = datasetlist->GetDataSet( analyzeArgs.GetStringKey("ds2") );
  if (check_type(ds2_,2)) return Analysis::ERR;

  mprintf("    DIVERGENCE: Between %s and %s\n", ds1_->Legend().c_str(),
          ds2_->Legend().c_str());

  return Analysis::OK;
}

// Analysis_Divergence::NormalizeSet()
std::vector<double> Analysis_Divergence::NormalizeSet(DataSet const& dsIn, unsigned int maxSize)
const {
  std::vector<double> setOut( maxSize );
  double sum = 0.0;
  DataSet_1D const& ds = static_cast<DataSet_1D const&>( dsIn );
  for (unsigned int i = 0; i < maxSize; i++)
    sum += ds.Dval(i);
  double norm = 1.0 / sum;
  for (unsigned int i = 0; i < maxSize; i++)
    setOut[i] = ds.Dval(i) * norm;
  return setOut;
}
  
// Analysis_Divergence::Analyze()
Analysis::RetType Analysis_Divergence::Analyze() {
  if (ds1_->Size() < 1 || ds2_->Size() < 1) {
    mprinterr("Error: One or both data sets empty (ds1=%zu, ds2=%zu)\n",
              ds1_->Size(), ds2_->Size());
    return Analysis::ERR;
  }
  size_t maxSize = ds1_->Size();
  if (maxSize != ds2_->Size()) {
    mprintf("Warning: Data set sizes do not match (ds1=%zu, ds2=%zu)\n",
              ds1_->Size(), ds2_->Size());
    maxSize = std::min( ds1_->Size(), ds2_->Size() );
    mprintf("Warning:  Only calculating divergence up to %zu points.\n", maxSize);
  }
  // Normalize sum over sets to 1.0
  std::vector<double> setP = NormalizeSet(*ds1_, maxSize);
  std::vector<double> setQ = NormalizeSet(*ds2_, maxSize);
  // Calculate divergence
  double divergence = 0.0;
  int nInvalid = 0;
  for (unsigned int i = 0; i < setP.size(); i++) {
    bool Pzero = (setP[i] <= std::numeric_limits<double>::min());
    bool Qzero = (setQ[i] <= std::numeric_limits<double>::min());
    if (!Pzero && !Qzero)
      divergence += (log( setP[i] / setQ[i] ) * setP[i]);
    else if ( Pzero != Qzero ) {
      // Only one value is zero, KL divergence not defined.
      nInvalid++;
    }
  }
  if (nInvalid > 0)
    mprintf("Warning:\tEncountered %i invalid points when calculating divergence.\n", nInvalid);
  mprintf("\tDivergence between %s and %s is %g\n", ds1_->Legend().c_str(),
          ds2_->Legend().c_str(), divergence);
  return Analysis::OK;
}
