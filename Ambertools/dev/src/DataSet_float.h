#ifndef INC_DATASET_FLOAT_H
#define INC_DATASET_FLOAT_H
#include <vector>
#include "DataSet_1D.h"
// Class: DataSet_float
/// Hold an array of float values.
class DataSet_float : public DataSet_1D {
  public:
    DataSet_float() : DataSet_1D(FLOAT, 8, 3) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_float();}
    float& operator[](size_t idx)        { return Data_[idx];         }
    float  operator[](size_t idx)  const { return Data_[idx];         }
    void AddElement(float f)             { Data_.push_back( f );      }
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)           { Data_.resize(sizeIn, 0.0); }
    // ----- DataSet functions -------------------
    size_t Size()                  const { return Data_.size();       }
    int Sync();
    void Info()                    const { return;                    }
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* );
    double Dval(size_t idx)        const { return (double)Data_[idx]; }
    double Xcrd(size_t idx)        const { return Dim(0).Coord(idx);  }
    void WriteBuffer(CpptrajFile&, size_t) const;
  private:
    std::vector<float> Data_;
};
#endif
