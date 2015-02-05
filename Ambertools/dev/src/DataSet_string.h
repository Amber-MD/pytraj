#ifndef INC_DATASET_STRING_H
#define INC_DATASET_STRING_H
#include <vector>
#include <string>
#include "DataSet_1D.h"
// Class: DataSet_string
/// Hold an array of string values.
class DataSet_string : public DataSet_1D {
  public:
    DataSet_string() : DataSet_1D(STRING, 1, 0) {}
    static DataSet* Alloc() { return (DataSet*)new DataSet_string();}
    std::string& operator[](size_t idx)  { return Data_[idx];         }
    void operator=(std::vector<std::string> const& rhs) { Data_ = rhs;}
    void AddElement(std::string const& s){ Data_.push_back( s );      }
    void Append(std::vector<std::string> const&);
    /// Make set size sizeIn, all values set to 0.0.
    void Resize(size_t sizeIn)           { Data_.resize(sizeIn, "");  }
    // ----- DataSet functions -------------------
    size_t Size()                  const { return Data_.size();       }
    int Sync();
    void Info()                    const { return;                    }
    // ----- DataSet_1D functions ----------------
    int Allocate1D(size_t);
    void Add( size_t, const void* );
    double Dval(size_t idx)        const { return 0.0;                }
    double Xcrd(size_t idx)        const { return Dim(0).Coord(idx);  }
    void WriteBuffer(CpptrajFile&, size_t) const;
  private:
    std::vector<std::string> Data_;
};
#endif
