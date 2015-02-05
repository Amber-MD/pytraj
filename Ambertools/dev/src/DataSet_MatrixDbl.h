#ifndef INC_DATASET_MATRIXDBL_H
#define INC_DATASET_MATRIXDBL_H
#include <vector>
#include "DataSet_2D.h"
#include "Matrix.h"
/// Double-precision two-dimensional matrix.
/** This is the class used by Action_Matrix. */
class DataSet_MatrixDbl : public DataSet_2D {
  public:
    DataSet_MatrixDbl() : DataSet_2D(MATRIX_DBL, 12, 4) {}
    double& operator[](size_t idx)             { return mat_[idx];          }
    static DataSet* Alloc() { return (DataSet*)new DataSet_MatrixDbl();     }
    // ----- DataSet functions -------------------
    size_t Size()                        const { return mat_.size();        }
    int Sync()                                 { return 1;                  }
    void Info()                          const { return;                    }
    // ----- DataSet_2D functions ----------------
    int Allocate2D(size_t x,size_t y)          { kind_=FULL; return mat_.resize(x,y); }
    int AllocateHalf(size_t x)                 { kind_=HALF; return mat_.resize(x,0); }
    int AllocateTriangle(size_t x)             { kind_=TRI;  return mat_.resize(0,x); }
    void Write2D(CpptrajFile&, int, int) const;
    double GetElement(size_t x,size_t y) const { return mat_.element(x,y);  }
    size_t Nrows()                       const { return mat_.Nrows();       }
    size_t Ncols()                       const { return mat_.Ncols();       }
    double* MatrixArray()                const;
    MatrixKind Kind()                    const { return kind_;              }
    MatrixType Type()                    const { return type_;              }
    // -------------------------------------------
    double& Element(size_t x, size_t y)        { return mat_.element(x,y);  }
    int AddElement(double d)                   { return mat_.addElement(d); }
    void SetElement(size_t x,size_t y,double d){ mat_.setElement(x,y,d);    }
    /// Type definition of iterator over matrix elements.
    typedef Matrix<double>::iterator iterator;
    iterator begin()                           { return mat_.begin();       }
    iterator end()                             { return mat_.end();         }
    /// Diagonal vector/mass array type definition. 
    typedef std::vector<double> Darray;
    /// \return diagonal vector.
    Darray const& Vect()                 const { return vect_;              }
    Darray& V1()                               { return vect_;              }
    /// Allocate diagonal vector.
    void AllocateVector(size_t vsize)          { vect_.resize(vsize, 0.0);  }
    /// \return iterator to beginning of diagonal vector.
    Darray::iterator v1begin()                 { return vect_.begin();      }
    /// \return iterator to end of diagonal vector.
    Darray::iterator v1end()                   { return vect_.end();        }
    /// Set matrix type and kind.
    void SetTypeAndKind(MatrixType tIn, MatrixKind kIn) { type_ = tIn; kind_ = kIn; }
    /// Store masses associated with columns in matrix.
    void StoreMass(Darray const& mIn)          { mass_ = mIn;               }
    /// \return array of masses associated with columns in matrix.
    Darray const& Mass()                 const { return mass_;              }
  private:
    Matrix<double> mat_;       ///< Matrix elements.
    Darray vect_;              ///< Hold diagonal elements | avg coords
    Darray mass_;              ///< Hold masses, for MWCOVAR quasiharmonic analysis. 
    MatrixType type_;          ///< Matrix type.
    MatrixKind kind_;          ///> Matrix kind.
};
#endif
