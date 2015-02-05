#ifndef INC_DATASET_DOUBLE_H
#define INC_DATASET_DOUBLE_H
#include <vector>
#include "DataSet_1D.h"
// Class: DataSet_double
/// Hold an array of double values.
class DataSet_double : public DataSet_1D {
  public:
    DataSet_double() : DataSet_1D(DOUBLE, 12, 4), bound_(0.0), 
                       boundh_(0.0), rexp_(-1.0) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_double();}
    double& operator[](size_t idx)       { return Data_[idx];         }
    double  operator[](size_t idx) const { return Data_[idx];         }
    std::vector<double> const& Data() const { return Data_;           }
    void operator=(std::vector<double> const& rhs) { Data_ = rhs;     }
    void AddElement(double d)            { Data_.push_back( d );      }
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)           { Data_.resize(sizeIn, 0.0); }
    // ----- DataSet functions -------------------
    size_t Size()                  const { return Data_.size();       }
    int Sync();
    void Info()                    const { return;                    }
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* );
    double Dval(size_t idx)        const { return Data_[idx];         }
    double Xcrd(size_t idx)        const { return Dim(0).Coord(idx);  }
    void WriteBuffer(CpptrajFile&, size_t) const;
    // -------------------------------------------
    void Append(std::vector<double> const&);
    void Append(DataSet_double const&);
    typedef std::vector<double>::iterator iterator;
    iterator begin()                     { return Data_.begin();      }
    iterator end()                       { return Data_.end();        }
    // For Analysis_Statistics DISTANCE NOE
    void SetNOE(double b, double bh, double r) { bound_=b; boundh_=bh; rexp_=r;}
    double NOE_bound()  const { return bound_;  }
    double NOE_boundH() const { return boundh_; }
    double NOE_rexp()   const { return rexp_;   }
    /// For dihedral/multidihedral/pucker, offset/wrap torsions.
    void ShiftTorsions(double,double);
  private:
    std::vector<double> Data_;
    // For Analysis_Statistics DISTANCE NOE
    double bound_;
    double boundh_;
    double rexp_;
};
#endif
