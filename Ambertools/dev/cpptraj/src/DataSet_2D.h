#ifndef INC_DATASET_2D_H
#define INC_DATASET_2D_H
#include "DataSet.h"
#include "CpptrajFile.h"
/// Interface for 2D DataSets.
class DataSet_2D : public DataSet {
  public:
    /// Types of matrix calculated by Action_Matrix.
    enum MatrixType {
      NO_OP=0, DIST, COVAR, MWCOVAR, CORREL, DISTCOVAR, IDEA, IRED, DIHCOVAR, NMAT
    };
    /// Kind of matrix.
    enum MatrixKind { FULL = 0, HALF, TRI };
    DataSet_2D() {}
    DataSet_2D(DataSet::DataType tIn, int wIn, int pIn) : 
      DataSet(tIn, wIn, pIn, 2) {}
    /// Set up matrix for given # rows and columns.
    virtual int Allocate2D(size_t, size_t) = 0;
    virtual int AllocateHalf(size_t) = 0;
    virtual int AllocateTriangle(size_t) = 0;
    /// Write 2D data to file (2D)
    virtual void Write2D(CpptrajFile&,int,int) const = 0;
    /// \return Data from matrix at col/row 
    virtual double GetElement(size_t, size_t) const = 0;
    /// \return the number of rows.
    virtual size_t Nrows() const = 0;
    /// \return the number of columns.
    virtual size_t Ncols() const = 0;
    /// \return double array containing matrix elements.
    virtual double* MatrixArray() const = 0;
    /// \return the kind of matrix, full/half/triangle.
    virtual MatrixKind Kind() const = 0;
    /// \return the type of matrix
    virtual MatrixType Type() const = 0;
    // -------------------------------------------
    // TODO: Remove this. Only needed by DataSet_1D.h
    void Add(size_t,const void*) { }
    static const char* MatrixTypeString(MatrixType m  ) { return TokenArray[m].TypeString;   }
    static const char* MatrixOutputString(MatrixType m) { return TokenArray[m].OutputString; }
  private:
    struct MatrixToken {
      const char* TypeString;
      const char* OutputString;
    };
    static const MatrixToken TokenArray[];
};
#endif
