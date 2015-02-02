#ifndef INC_DATASET_H
#define INC_DATASET_H
#include <cstddef> // size_t
#include <vector>
#include "Dimension.h"
// Class: DataSet
/// Base class that all DataSet types will inherit.
/** DataSets are given certain attributes to make DataSet selection easier; 
  * these are name, index, and aspect. Name is typically associated with the
  * action that creates the DataSet, e.g. RMSD or distance. Index is used
  * when and action outputs subsets of data, e.g. with RMSD it is possible to 
  * output per-residue RMSD, where the DataSet index corresponds to the residue
  * number. Aspect is used to further subdivide output data type; e.g. with 
  * nucleic acid analysis each base pair (divided by index) has shear,
  * stagger etc calculated.
  */
class DataSet {
  public:
    typedef DataSet* (*AllocatorType)();
    /// Type of data stored in DataSet
    enum DataType {
      UNKNOWN_DATA=0, DOUBLE, FLOAT, INTEGER, STRING, MATRIX_DBL, MATRIX_FLT, 
      COORDS, VECTOR, MODES, GRID_FLT, REMLOG, XYMESH, TRAJ
    };
    /// Source of data stored in DataSet, used by Analysis_Statistics
    enum scalarMode {
      M_DISTANCE=0, M_ANGLE, M_TORSION, M_PUCKER, M_RMS, UNKNOWN_MODE 
    };
    /// Type of DataSet, used by Analysis_Statistics
    enum scalarType {
      ALPHA=0, BETA, GAMMA, DELTA, EPSILON, ZETA,  PUCKER, CHI, H1P,
      C2P,     PHI,  PSI,   PCHI,  OMEGA,   NOE,   UNDEFINED
    };

    DataSet();
    /// Set DataSet type, width, precision, and dimension.
    DataSet(DataType,int,int,int);
    DataSet(const DataSet&);
    DataSet& operator=(const DataSet&);
    virtual ~DataSet() {} // Destructor - virtual since this class is inherited

    // ----------===== Inheritable functions =====----------
    /// \return the number of data elements stored in the DataSet.
    virtual size_t Size() const = 0;
    /// Consolidate this DataSet across all threads (MPI only)
    virtual int Sync() = 0;
    /// Print DataSet information
    virtual void Info() const = 0;
    // -----------------------------------------------------
    // TODO: Remove this. Should only be in DataSet_1D.h
    virtual void Add( size_t, const void* ) = 0;
    // ----------===== Public functions =====---------------
    /// Set output width.
    void SetWidth(int);
    /// Set output width and precision
    void SetPrecision(int,int);
    /// Set up DataSet with given name, index, and aspect.
    int SetupSet(std::string const&,int,std::string const&);
    /// Set DataSet legend.
    void SetLegend( std::string const& lIn ) { legend_ = lIn;     }
    /// Set scalar mode
    void SetScalar( scalarMode mIn )         { scalarmode_ = mIn; }
    /// Set specified DataSet dimension.
    void SetDim(Dimension::DimIdxType i, Dimension const& d) { dim_[(int)i]=d; }
    /// Set scalar mode and type
    inline void SetScalar( scalarMode, scalarType );
    /// Used to set the data and header format strings 
    int SetDataSetFormat(bool);
    /// Check if name and/or index and aspect match this DataSet.
    bool Matches(std::string const&, int, std::string const&) const;
    /// Write scalar mode/type description
    void ScalarDescription() const;
    // -----------------------------------------------------
    // ---===== Functions that return private vars =====----
    /// True if DataSet is empty. 
    bool Empty()                const { return (Size() == 0);      }
    /// DataSet output label.
    std::string const& Legend() const { return legend_;            }
    // TODO: Put all data set metadata into one class.
    /// \return DataSet base name.
    std::string const& Name()   const { return name_;              }
    /// \return DataSet index.
    int Idx()                   const { return idx_;               }
    /// \return DataSet aspect.
    std::string const& Aspect() const { return aspect_;            }
    /// \return Total output width of a data element in characters.
    int ColumnWidth()           const { return colwidth_;          }
    /// \return DataSet type.
    DataType Type()             const { return dType_;             }
    /// \return scalar mode
    scalarMode ScalarMode()     const { return scalarmode_;        }
    /// \return scalar type
    scalarType ScalarType()     const { return scalartype_;        }
    /// \return number of dimensions.
    size_t Ndim()               const { return dim_.size();        }
    /// \return true if set is a coordinate set type.
    bool IsCoordSet()           const { return (dType_ == COORDS || dType_ == TRAJ); }
    /// \return specified DataSet dimension.
    Dimension& Dim(Dimension::DimIdxType i) { return dim_[(int)i]; }
    Dimension&       Dim(int i)             { return dim_[i];      }
    Dimension const& Dim(int i)       const { return dim_[i];      }
    /// Comparison for sorting, name/aspect/idx
    inline bool operator<(const DataSet&) const;
    /// Used to sort DataSet pointers (DataSetList, Array1D).
    struct DS_PtrCmp {
      inline bool operator()(DataSet const* first, DataSet const* second) const {
        return *first < *second;
      }
    };
    const char* DataFormat()    const { return data_format_;       }
    static const char* Smodes[];
    static const char* Stypes[];
    static const scalarMode TypeModes[];
    /// \return scalarMode that matches input keyword.
    static scalarMode ModeFromKeyword(std::string const&);
    /// \return scalarType that matches keyword; check that mode is valid if specified.
    static scalarType TypeFromKeyword(std::string const&, scalarMode&);
    static scalarType TypeFromKeyword(std::string const&, scalarMode const&);
  protected:
    /// Width of numbers in output elements.
    int Width()                 const { return width_;             }
    const char* data_format_; ///< Used to avoid constant calls to format_.c_str().
  private:
    /// Type to hold coordinate info for each dimension in DataSet.
    typedef std::vector<Dimension> DimArray;
    std::string name_;        ///< Name of the DataSet
    int idx_;                 ///< DataSet index
    std::string aspect_;      ///< DataSet aspect.
    std::string legend_;      ///< DataSet legend.
    DataType dType_;          ///< The DataSet type
    DimArray dim_;            ///< Holds info for each dimension in the DataSet.
    int colwidth_;            ///< The total output width of a data element.
    int width_;               ///< The output width of numbers in a data element.
    int precision_;           ///< The output precision of numbers in a data element.
    bool leftAlign_;          ///< If true output will be left-aligned (no leading space).
    std::string format_;      ///< Output printf format string for data.
    scalarMode scalarmode_;   ///< Source of data in DataSet.
    scalarType scalartype_;   ///< Specific type of data in DataSet (if any).
};
// ---------- INLINE FUNCTIONS -------------------------------------------------
bool DataSet::operator<(const DataSet& rhs) const {
  if ( name_ == rhs.name_ ) {
    if ( aspect_ == rhs.aspect_ ) {
      return ( idx_ < rhs.idx_ );
    } else {
      return ( aspect_ < rhs.aspect_ );
    }
  } else {
    return ( name_ < rhs.name_ );
  }
}

void DataSet::SetScalar( scalarMode modeIn, scalarType typeIn ) {
  scalarmode_ = modeIn;
  scalartype_ = typeIn;
}
#endif 
