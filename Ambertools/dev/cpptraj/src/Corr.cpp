#include "Corr.h"
// CorrF_Direct::Allocate()
void CorrF_Direct::Allocate(int stepsIn) {
  nsteps_ = stepsIn;
  table_.resize(2*nsteps_, 0.0);
}

// CorrF_Direct::AutoCorr()
void CorrF_Direct::AutoCorr(ComplexArray& data1) {
  int ndata2 = data1.size();
  for (int i = 0; i < ndata2; i++) {
    double dsum = 0.0;
    for (int j = i; j < ndata2; j++) {
      int ind1 = 2 * j;
      int ind2 = 2 * (j-i);
      dsum += data1[ind1] * data1[ind2] + data1[ind1+1] * data1[ind2+1];
    }
    if (i < nsteps_) {
      int ind1 = 2 * i;
      table_[ind1  ] = dsum;
      table_[ind1+1] = 0.0;
    } else
      break;
  }
  std::copy(table_.begin(), table_.end(), data1.CAptr());
}

// CorrF_Direct::CrossCorr()
void CorrF_Direct::CrossCorr(ComplexArray& data1, ComplexArray const& data2) {
  if (data2.size() < data1.size()) return;
  int ndata2 = data1.size();
  for (int i = 0; i < ndata2; i++) {
    double dsum = 0.0;
    double dsumi = 0.0;
    for (int j = i; j < ndata2; j++) {
      int ind1 = 2 * j;
      int ind2 = 2 * (j-i);
      dsum  += data2[ind1] * data1[ind2  ] + data2[ind1+1] * data1[ind2+1];
      dsumi += data2[ind1] * data1[ind2+1] - data2[ind1+1] * data1[ind2  ];
    }
    if(i < nsteps_) {
      int ind1 = 2 * i;
      table_[ind1  ] = dsum;
      table_[ind1+1] = dsumi;
    } else
      break;
  }
  std::copy(table_.begin(), table_.end(), data1.CAptr());
}

// -----------------------------------------------------------------------------
// CorrF_FFT::Allocate()
void CorrF_FFT::Allocate(int stepsIn) {
  pubfft_.SetupFFT_NextPowerOf2( stepsIn );
}

// CorrF_FFT::AutoCorr()
void CorrF_FFT::AutoCorr(ComplexArray& data1) {
  pubfft_.Forward( data1 );
  // Calculate square modulus of F(data1)
  data1.SquareModulus();
  // Inverse FFT
  pubfft_.Back( data1 );
  // Normalize with fft_size (since not done in inverse FFT routine)
  data1.Normalize( 1.0 / (double)pubfft_.size() );
}

// CorrF_FFT::CrossCorr()
void CorrF_FFT::CrossCorr(ComplexArray& data1, ComplexArray& data2) {
  // Cross-correlation
  pubfft_.Forward( data1 );
  pubfft_.Forward( data2 );
  // Calculate [data1]* x [data2] where * denotes complex conjugate.
  data1.ComplexConjTimes(data2);
  // Inverse FFT
  pubfft_.Back( data1 );
  // Normalize with fft_size (since not done in inverse FFT routine)
  data1.Normalize( 1.0 / (double)pubfft_.size() );
}
