// Collection of routines to perform math on 1D datasets.
#include <cmath> // sqrt, fabs
#include "DataSet_1D.h"
#include "Corr.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD, RADDEG

/// Return true if set is an atomic type (i.e. int, double, float).
bool DataSet_1D::GoodCalcType(DataSet_1D const& ds) {
  if (ds.Type() == DataSet::DOUBLE  ||
      ds.Type() == DataSet::FLOAT   ||
      ds.Type() == DataSet::INTEGER ||
      ds.Type() == DataSet::XYMESH    )
    return true;
  mprinterr("Error: DataSet %s is not a valid type for this calc.\n",
            ds.Name().c_str());
  return false;
}

/** Calculate the average over values in this set (and optionally the
  * standard deviation).
  */
double DataSet_1D::Avg(double* stdev) const {
  // Check # values
  int numvalues = Size();
  if ( numvalues < 1 ) {
    if (stdev!=0) *stdev = 0.0;
    return 0.0;
  }
  double avg = 0;
  // Check if this set is a good type
  if ( GoodCalcType(*this) ) {
    if (IsTorsionArray()) {
      // Cyclic torsion average
      double sumy = 0.0;
      double sumx = 0.0;
      for ( int i = 0; i < numvalues; ++i ) {
        double theta = Dval( i ) * Constants::DEGRAD;
        sumy += sin( theta );
        sumx += cos( theta );
      }
      avg = atan2(sumy, sumx) * Constants::RADDEG;
      if (stdev==0) return avg;
      // Torsion Stdev
      sumy = 0;
      for ( int i = 0; i < numvalues; ++i) {
        double diff = fabs(avg - Dval( i ));
        if (diff > 180.0)
          diff = 360.0 - diff;
        diff *= diff;
        sumy += diff;
      }
      sumy /= (double)numvalues;
      *stdev = sqrt(sumy);
    } else {
      // Non-cyclic, normal average
      double sum = 0;
      for ( int i = 0; i < numvalues; ++i )
        sum += Dval( i );
      avg = sum / (double)numvalues;
      if (stdev==0) return avg;
      // Stdev
      sum = 0;
      for ( int i = 0; i < numvalues; ++i ) {
        double diff = avg - Dval( i );
        diff *= diff;
        sum += diff;
      }
      sum /= (double)numvalues;
      *stdev = sqrt(sum);
    }
  }
  return avg;
}

/** Return the minimum value in the dataset.  */
double DataSet_1D::Min() const {
  // Check # values
  if (Size()==0) return 0;
  double min = 0;
  // Check if this set is a good type
  if ( GoodCalcType(*this) ) {
    min = Dval( 0 );
    for (size_t i = 1; i < Size(); ++i) {
      double val = Dval( i );
      if (val < min) min = val;
    }
  }
  return min;
}

/** Return the maximum value in the dataset.  */
double DataSet_1D::Max() const {
  // Check # values
  if ( Size() == 0 ) return 0;
  double max = 0;
  // Check if this set is a good type
  if ( GoodCalcType(*this) ) {
    max = Dval( 0 );
    for (size_t i = 1; i < Size(); ++i) {
      double val = Dval( i );
      if (val > max) max = val;
    }
  }
  return max;
}

static inline double PeriodicDiff(double v1, double v2) {
  double diff = v1 - v2;
  if (diff > 180.0)
    diff = 360.0 - diff;
  else if (diff < -180.0)
    diff = 360.0 + diff;
  return diff;
}

/** Calculate time correlation between two DataSets.
  * \D2 DataSet to calculate correlation to.
  * \Ct DataSet to store time correlation fn, must be DOUBLE.
  * \lagmaxIn Max lag to calculate corr. -1 means use size of dataset.
  * \calccovar If true calculate covariance (devation from avg).
  * \return 0 on success, 1 on error.
  */
int DataSet_1D::CrossCorr( DataSet_1D const& D2, DataSet_1D& Ct,
                           int lagmaxIn, bool calccovar, bool usefft ) const
{
  int lagmax;
  double ct;
  // Check if D1 and D2 are valid types
  if ( !GoodCalcType(*this) ) return 1;
  if ( !GoodCalcType(D2) ) return 1;
  // Check that D1 and D2 have same # data points.
  // TODO: size_t
  int Nelements = (int)Size();
  if (Nelements != (int)D2.Size()) {
    mprinterr("Error: CrossCorr: # elements in dataset %s (%i) not equal to\n",
              Legend().c_str(), Nelements);
    mprinterr("Error:            # elements in dataset %s (%u)\n",
              D2.Legend().c_str(), D2.Size());
    return 1;
  }
  if (Nelements < 2) {
    mprinterr("Error: CrossCorr: # elements is less than 2 (%i)\n", Nelements);
    return 1;
  }
  // Check return dataset type
  if ( Ct.Type() != DataSet::DOUBLE ) {
    mprinterr("Internal Error: CrossCorr: Ct must be of type DataSet::DOUBLE.\n");
    return 1;
  }
  // Check if lagmaxIn makes sense. Set default lag to be Nelements 
  // if not specified.
  if (lagmaxIn == -1)
    lagmax = Nelements;
  else if (lagmaxIn > Nelements) {
    mprintf("Warning: CrossCorr [%s][%s]: max lag (%i) > Nelements (%i), setting to Nelements.\n",
            Legend().c_str(), D2.Legend().c_str(), lagmaxIn, Nelements);
    lagmax = Nelements;
  } else
    lagmax = lagmaxIn;
  // If calculating covariance calculate averages
  double avg1 = 0;
  double avg2 = 0;
  if ( calccovar ) {
    avg1 = this->Avg();
    avg2 = D2.Avg();
  }
  // Calculate correlation
  double norm = 1.0;
  if ( usefft ) {
    // Calc using FFT
    CorrF_FFT pubfft1(Nelements);
    ComplexArray data1 = pubfft1.Array();
    data1.PadWithZero(Nelements);
    if (IsTorsionArray()) {
      for (int i = 0; i < Nelements; ++i)
        data1[i*2] = PeriodicDiff(avg1, Dval( i ));
    } else {
      for (int i = 0; i < Nelements; ++i)
        data1[i*2] = Dval(i) - avg1;
    }
    if (&D2 == this)
      pubfft1.AutoCorr(data1);
    else {
      // Populate second dataset if different
      ComplexArray data2 = pubfft1.Array();
      data2.PadWithZero(Nelements);
      if (D2.IsTorsionArray()) {
        for (int i = 0; i < Nelements; ++i)
          data2[i*2] = PeriodicDiff(avg2, D2.Dval( i ));
      } else {
        for (int i = 0; i < Nelements; ++i)
          data2[i*2] = D2.Dval(i) - avg2;
      }
      pubfft1.CrossCorr(data1, data2);
    }
    // Put real components of data1 in output DataSet
    norm = 1.0 / fabs( data1[0] );
    for (int i = 0; i < lagmax; ++i) {
      ct = data1[i*2] * norm;
      Ct.Add(i, &ct);
    }
  } else {
    // Direct calc
    double diff1, diff2;
    for (int lag = 0; lag < lagmax; ++lag) {
      ct = 0;
      int jmax = Nelements - lag;
      for (int j = 0; j < jmax; ++j) {
        if (IsTorsionArray())
          diff1 = PeriodicDiff(Dval(j), avg1);
        else
          diff1 = Dval(j) - avg1;
        if (D2.IsTorsionArray())
          diff2 = PeriodicDiff(D2.Dval(j+lag), avg2);
        else   
          diff2 = D2.Dval(j+lag) - avg2;
        ct += (diff1 * diff2);
      }
      if (lag == 0) {
        if (ct != 0)
          norm = fabs( ct );
      }
      ct /= norm;
      Ct.Add(lag, &ct);
    }
  }
  return 0;
}

/** Calculate Pearson product-moment correlation between DataSets.
  * \D2 DataSet to caclulate correlation to.
  * \return Pearson product-moment correlation coefficient.
  */
double DataSet_1D::CorrCoeff( DataSet_1D const& D2 ) const {
  // Check if D1 and D2 are valid types
  if ( !GoodCalcType(*this) ) return 0;
  if ( !GoodCalcType(D2) ) return 0;
  // Check that D1 and D2 have same # data points.
  // TODO: size_t
  int Nelements = (int)Size();
  if (Nelements != (int)D2.Size()) {
    mprinterr("Error: Corr: # elements in dataset %s (%i) not equal to\n",
              Legend().c_str(), Nelements);
    mprinterr("Error:       # elements in dataset %s (%u)\n",
              D2.Legend().c_str(), D2.Size());
    return 0;
  }
  // Calculate averages
  double avg1 = this->Avg();
  double avg2 = D2.Avg();
  // Calculate average deviations. 
  double sumdiff1_2 = 0.0;
  double sumdiff2_2 = 0.0;
  double corr_coeff = 0.0;
  //mprinterr("DATASETS %s and %s\n", c_str(), D2.c_str());
  for (int i = 0; i < Nelements; i++) {
    double diff1 = Dval(i) - avg1;
    double diff2 = D2.Dval(i) - avg2;
    sumdiff1_2 += (diff1 * diff1);
    sumdiff2_2 += (diff2 * diff2);
    corr_coeff += (diff1 * diff2);
  }
  if (sumdiff1_2 == 0.0 || sumdiff2_2 == 0.0) {
    mprintf("Warning: Corr: %s to %s, Normalization is 0\n",
            Legend().c_str(),  D2.Legend().c_str());
    return 0;
  }
  // Correlation coefficient
  corr_coeff /= ( sqrt( sumdiff1_2 ) * sqrt( sumdiff2_2 ) );
  //mprintf("    CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
  //        D1_->c_str(), D2_->c_str(), corr_coeff );
  return corr_coeff;
}
