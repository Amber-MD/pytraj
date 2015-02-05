#ifndef INC_COMPLEXARRAY_H
#define INC_COMPLEXARRAY_H
#include "ArrayIterator.h"
/// Array that will hold complex numbers.
/** Implementation does not use STL vector so as to easily interface with
  * fortran routines for e.g. calculating FFT; current standard does not
  * directly allow access interal STL vector representation (i.e. double*).
  */
class ComplexArray {
  public:
    ComplexArray() : data_(0), ndata_(0), ncomplex_(0) {}
    ComplexArray(int);
    ComplexArray(ComplexArray const&);
    ComplexArray& operator=(ComplexArray const&);
    ~ComplexArray();
    void Allocate(int);
    /// Copy data from dataIn but do not reallocate.
    void Assign(ComplexArray const&);
    /// Pad array with zero from given position to end.
    void PadWithZero(int);
    /// Multiply all entries in data by input 
    void Normalize(double);
    void SquareModulus();
    /// Calculate [this]* x [rhs] where * denotes complex conjugate.
    void ComplexConjTimes(ComplexArray const&);
    double* CAptr()  { return data_;     }
    int size() const { return ncomplex_; }
    double& operator[](int idx)             { return data_[idx]; }
    double const& operator[](int idx) const { return data_[idx]; }
    // iterators
    typedef ArrayIterator<double> iterator;
    const iterator begin() const { return data_;          }
    const iterator end()   const { return data_ + ndata_; }
  private:
    double* data_;
    int ndata_;    ///< Actual size of data
    int ncomplex_; ///< # of complex numbers in array (ndata / 2)
};
#endif
