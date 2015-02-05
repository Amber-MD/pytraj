#ifndef INC_ARRAY1D_H
#define INC_ARRAY1D_H
#include "DataSetList.h"
#include "DataSet_1D.h"
#include "ArgList.h"
/// Hold 1D DataSets
class Array1D {
  public:
    Array1D() {}
    Array1D(const Array1D&);
    Array1D& operator=(const Array1D&);
    Array1D(DataSetList const&);
    /// \return size of largest data set in the array.
    size_t DetermineMax() const;
    int push_back( DataSet_1D* const& );
    DataSet_1D* const& operator[](int idx) const { return array_[idx];       }
    DataSet_1D* operator[](int idx)              { return array_[idx];       }
    bool empty()                           const { return array_.empty();    }
    typedef std::vector<DataSet_1D*>::const_iterator const_iterator;
    const_iterator begin()                 const { return array_.begin();    }
    const_iterator end()                   const { return array_.end();      }
    size_t size()                          const { return array_.size();     }
    void clear()                                 { array_.clear();           }
    void SortArray1D();
    int AddDataSets(DataSetList const&);
    int AddTorsionSets(DataSetList const&);
    int AddSetsFromArgs(ArgList const&, DataSetList const&);
    int CheckXDimension() const;
  private:
    std::vector<DataSet_1D*> array_;
};
#endif
