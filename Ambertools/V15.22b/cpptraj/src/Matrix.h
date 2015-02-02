#ifndef INC_MATRIX_H
#define INC_MATRIX_H
#include "ArrayIterator.h"
/// Two-dimensional matrix template.
template <class T> class Matrix {
// TODO: Type may not be necessary here if in DataSet_2D
    enum MType { FULL = 0, HALF, TRIANGLE }; 
  public:
    Matrix() :
      elements_(0), ncols_(0L), nrows_(0L), nelements_(0L),
      maxElements_(0L), currentElement_(0L), type_(FULL),
      calcIndex(calcFullIndex) {}
    ~Matrix() { if (elements_!=0) delete[] elements_; }
    Matrix( const Matrix& );
    Matrix& operator=( const Matrix& );
    T& operator[](size_t idx)                  { return elements_[idx];  }
    const T& operator[](size_t idx)      const { return elements_[idx];  }
    /// \return total number of elements in the matrix.
    size_t size()                        const { return nelements_;      }
    /// \return current matrix type.
    //MType Type()                         const { return type_;           }
    /// Set up matrix for given number of cols and rows.
    int resize(size_t,size_t);
    /// \return element at specified col and row.
    const T& element(int,int) const;
    T&       element(int,int);
    /// \return number of rows (Y).
    size_t Nrows()     const { return nrows_;     }
    /// \return number of cols (X).
    size_t Ncols()     const { return ncols_;     }
    /// Add an element to the matrix in order.
    int addElement( const T& );
    /// Set element at col and row.
    void setElement(int,int, const T&);
    /// \return pointer to internal array of elements.
    T const* Ptr()     const { return elements_;  }
    T* Ptr()                 { return elements_;  }
    /// Convert X and Y to index. 
    size_t CalcIndex(int x, int y) const { return calcIndex(ncols_, x, y); }
    /// Iterator over matrix elements
    typedef ArrayIterator<T> iterator;
    /// Iterator to beginning of matrix elements.
    iterator begin() { return elements_;              }
    iterator begin() const { return elements_; }
    /// Iterator to end of matrix elements.
    iterator end()   { return elements_ + nelements_; }
    iterator end()   const { return elements_ + nelements_; }
    /// Return memory used by matrix in bytes.
    size_t DataSize() const {
      return (nelements_*sizeof(T)) + sizeof(T) +
             (5 * sizeof(size_t) + sizeof(MType) +
             sizeof(long int(*)()));
    }
  private:
    T* elements_;           ///< Array of elements
    T diagElt_;             ///< For TRIANGLE, the value of the diagonal.
    size_t ncols_;          ///< Number of columns (X)
    size_t nrows_;          ///< Number of rows (Y)
    size_t nelements_;      ///< Total number of elements.
    size_t maxElements_;    ///< Max number of elements currently allocated for.
    size_t currentElement_; ///< Current element (for AddElement())
    MType type_;            ///< Current matrix type.
    /// Pointer to index calculator for current matrix type
    long int (*calcIndex)(size_t,int,int);
    /// Full 2D matrix. 
    static long int calcFullIndex(size_t nX,int x,int y) {
      return (long int)( (y*(int)nX)+x );
    }
    /// Upper-half matrix + diagonal.
    static long int calcHalfIndex(size_t nX, int xIn, int yIn) {
      int i, j;
      if (yIn > xIn) {
        i = xIn;
        j = yIn;
      } else {
        i = yIn;
        j = xIn;
      }
      return (long int)(i * (int)nX - (i * (i-1) / 2) + (j - i));
    }
    /// Upper-half matrix - diagonal.
    static long int calcTriIndex(size_t nX, int xIn, int yIn) {
      int i, j;
      if (yIn > xIn) {
        i = xIn;
        j = yIn;
      } else if (xIn > yIn) {
        i = yIn;
        j = xIn;
      } else // iIn == jIn, triangle matrix diagonal is indicated by -1 
        return -1L;
      int i1 = i + 1;
      return ( ( ((int)nX * i) - ((i1 * i) / 2) ) + j - i1 );
    }
};
// COPY CONSTRUCTOR
template<class T> Matrix<T>::Matrix(const Matrix& rhs) :
  elements_(0),
  diagElt_( rhs.diagElt_ ),
  ncols_( rhs.ncols_ ),
  nrows_( rhs.nrows_ ),
  nelements_( rhs.nelements_ ),
  maxElements_( rhs.maxElements_ ),
  currentElement_( rhs.currentElement_ ),
  type_( rhs.type_ ),
  calcIndex( rhs.calcIndex )
{
  if (maxElements_ > 0L) {
    elements_ = new T[ maxElements_ ];
    std::copy( rhs.elements_, rhs.elements_ + nelements_, elements_ );
  }
}
// ASSIGNMENT
template<class T> Matrix<T>& Matrix<T>::operator=(const Matrix& rhs) {
  if (this == &rhs) return *this;
  if (elements_!=0) {
    delete[] elements_;
    elements_ = 0;
  }
  ncols_ = rhs.ncols_;
  nrows_ = rhs.nrows_;
  nelements_ = rhs.nelements_;
  maxElements_ = rhs.maxElements_;
  diagElt_ = rhs.diagElt_;
  currentElement_ = rhs.currentElement_;
  type_  = rhs.type_;
  calcIndex = rhs.calcIndex;
  if (maxElements_ > 0L) {
    elements_ = new T[ maxElements_ ];
    std::copy( rhs.elements_, rhs.elements_ + nelements_, elements_ );
  }
  return *this;
}
// Matrix::resize()
/** \param nX number of columns.
  * \param nY number of rows.
  * Three types of allocation:
  *   - cols > 0, rows > 0: Full matrix
  *   - cols > 0, rows == 0: Half matrix
  *   - cols == 0, rows > 0: Triangle matrix
  */
template<class T> int Matrix<T>::resize(size_t nX, size_t nY) {
  diagElt_ = T(); // Diagonal element default to zero.
  if (nX > 0L && nY > 0L) { // FULL
    ncols_ = nX;
    nrows_ = nY;
    nelements_ = ncols_ * nrows_;
    calcIndex = calcFullIndex;
    type_ = FULL;
  } else if (nX > 0L && nY == 0L) { // HALF
    ncols_ = nX;
    nrows_ = nX;
    nelements_ = ncols_ * (ncols_ + 1L) / 2L;
    calcIndex = calcHalfIndex;
    type_ = HALF;
  } else if (nX == 0L && nY > 0L) { // TRIANGLE
    ncols_ = nY;
    nrows_ = nY;
    nelements_ = (ncols_ * (ncols_ - 1L) / 2L);
    calcIndex = calcTriIndex;
    type_ = TRIANGLE;
  } else { // Both Zero, EMPTY
    ncols_ = 0L;
    nrows_ = 0L;
    nelements_ = 0L;
    return 1;
  }
  currentElement_ = 0L;
  if (nelements_ > 0L) {
    if (nelements_ > maxElements_) {
      if (elements_ != 0) delete[] elements_;
      elements_ = new T[ nelements_ ];
      maxElements_ = nelements_;
    }
    std::fill(elements_, elements_ + nelements_, T());
  }
  return 0;
}
// Matrix::addElement()
/** Add the input T to the element array and increment currentElement.
  * \return 1 on success, 0 if no more elements can be added.
  */
template<class T> int Matrix<T>::addElement(const T& elementIn) {
  if (currentElement_>=nelements_) return 0;
  elements_[currentElement_] = elementIn;
  ++currentElement_;
  return 1;
}
// Matrix::setElement()
template<class T> void Matrix<T>::setElement(int xIn, int yIn, const T& eltIn) {
  long int idx = calcIndex(ncols_, xIn, yIn);
  elements_[idx] = eltIn;
}
// Matrix::element()
template<class T> const T& Matrix<T>::element(int xIn, int yIn) const {
  long int idx = calcIndex(ncols_, xIn, yIn);
  if (idx < 0) return diagElt_; // In case of xIn == yIn for TRIANGLE
  return elements_[idx];
}
template<class T> T& Matrix<T>::element(int xIn, int yIn) {
  long int idx = calcIndex(ncols_, xIn, yIn);
  if (idx < 0) return diagElt_; // In case of xIn == yIn for TRIANGLE
  return elements_[idx];
}
#endif
