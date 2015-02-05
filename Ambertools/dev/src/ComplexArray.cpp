#include <algorithm> // fill, copy
#include "ComplexArray.h"

// CONSTRUCTOR
ComplexArray::ComplexArray(int i) : ndata_(i*2), ncomplex_(i)
{ 
  if (ndata_ > 0) { 
    data_ = new double[ ndata_ ];
    std::fill(data_, data_ + ndata_, 0);
  } else
    data_ = 0;
}

// COPY CONSTRUCTOR
ComplexArray::ComplexArray(ComplexArray const& rhs) :
  ndata_(rhs.ndata_),
  ncomplex_(rhs.ncomplex_)
{
  if (ndata_ > 0) { 
    data_ = new double[ ndata_ ];
    std::copy(rhs.data_, rhs.data_ + ndata_, data_);
  } else
    data_ = 0;
}

// ASSIGNMENT
ComplexArray& ComplexArray::operator=(ComplexArray const& rhs) {
  if (&rhs == this) return *this;
  if (data_ != 0) delete[] data_;
  ndata_ = rhs.ndata_;
  ncomplex_ = rhs.ncomplex_;
  if (ndata_ > 0) {
    data_ = new double[ ndata_ ];
    std::copy(rhs.data_, rhs.data_ + ndata_, data_);
  } else
    data_ = 0;
  return *this;
}

// DESTRUCTOR
ComplexArray::~ComplexArray() { if (data_ != 0) delete[] data_; }

// ComplexArray::Assign()
void ComplexArray::Assign(ComplexArray const& dataIn) {
  // TODO: Check that ndata >= dataIn.ndata
  std::copy( dataIn.data_, dataIn.data_ + dataIn.ndata_, data_ );
}

// ComplexArray::Allocate()
void ComplexArray::Allocate(int i) {
  ndata_ = i * 2;
  ncomplex_ = i;
  if (data_!=0) delete[] data_;
  if (ndata_ > 0) {
    data_ = new double[ ndata_ ];
    std::fill(data_, data_ + ndata_, 0);
  } else
    data_ = 0;
}

// ComplexArray::PadWithZero()
void ComplexArray::PadWithZero(int start) {
  int startIdx = 2 * start;
  std::fill(data_ + startIdx, data_ + ndata_, 0);
}

// ComplexArray::Normalize()
void ComplexArray::Normalize(double norm) {
  for (int i = 0; i < ndata_; ++i)
    data_[i] *= norm;
}

// ComplexArray::SquareModulus()
void ComplexArray::SquareModulus() {
  for (int i = 0; i < ndata_; i+=2) {
    data_[i  ] = data_[i  ] * data_[i  ] + data_[i+1] * data_[i+1];
    data_[i+1] = 0.0;
  }
}

// ComplexArray::ComplexConjTimes()
// TODO: Check size of rhs
void ComplexArray::ComplexConjTimes(ComplexArray const& rhs) {
  for (int i = 0; i < ndata_; i+=2) {
    double dtmp = data_[i  ] * rhs.data_[i  ] + data_[i+1] * rhs.data_[i+1];
    data_[i+1]  = data_[i  ] * rhs.data_[i+1] - data_[i+1] * rhs.data_[i  ];
    data_[i  ]  = dtmp;
  }
}
