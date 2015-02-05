#ifndef INC_DATASETLIST_H
#define INC_DATASETLIST_H
#include <vector>
#include "DataSet.h"
#include "ArgList.h" // GetReferenceFrame
#include "ReferenceFrame.h" // GetReferenceFrame
// Class: DataSetList
/// Hold list of data sets.
/** Main class for handling datasets. All dataset types can be allocated 
  * by DataSetList. DataSets are added to the list by various actions. 
  * There is a master DataSetList in Cpptraj that will hold data indexed 
  * by frame for use in analysis. All Data that may need to be available
  * for analysis should be in the master DataSetList.
  * This class is also used by DataFile to hold DataSets to be printed out.
  * This is done primarily to make use of the Get() functionality of 
  * DataSetList. In DataFile only copies of the sets are held, so they
  * are not freed by the destructor.
  */
class DataSetList {
  public:
    DataSetList();
    ~DataSetList();
    void Clear();
    DataSetList& operator+=(DataSetList const&);
    /// \return Description for given set type.
    static const char* SetString(DataSet::DataType);
    /// DataSetList default iterator
    typedef std::vector<DataSet*>::const_iterator const_iterator;
    /// Iterator to beginning of dataset list
    const_iterator begin() const { return DataList_.begin(); }
    /// Iterator to end of dataset list
    const_iterator end()   const { return DataList_.end();   }
    /// True if no DataSets in list.
    bool empty()           const { return DataList_.empty(); }
    /// \return number of datasets in the list 
    size_t size()          const { return DataList_.size();  }
    /// \return Ensemble number; -1 if not an ensemble
    int EnsembleNum()      const { return ensembleNum_;      }
    /// \return True if Actions have indicated DataSets will be generated.
    bool DataSetsPending() const { return dataSetsPending_;  }
    /// Remove set from list - used in DataFile
    void RemoveSet( const_iterator );
    /// Remove set from the list.
    void RemoveSet( DataSet* );
    /// \return DataSet at didx.
    DataSet* operator[](int didx) { return DataList_[didx]; } // FIXME: No bounds check
    /// Set DataSetList and underlying DataSet debug level
    void SetDebug(int);
    /// Set current ensemble number.
    void SetEnsembleNum(int i)   { ensembleNum_ = i;        }
    /// Set DataSets pending status.
    void SetDataSetsPending(bool b) { dataSetsPending_ = b; }
    /// Allocate 1D DataSet memory based on current max# expected frames.
    void AllocateSets(long int);
    /// Set width.precision of all DataSets in the list.
    void SetPrecisionOfDataSets(std::string const&, int, int);
    /// Get DataSet with specified name, index, and aspect.
    DataSet* GetSet(std::string const&, int, std::string const&) const;
    /// Get DataSet matching specified argument.
    DataSet* GetDataSet( std::string const& ) const;
    /// Get DataSet matching specified argument, no warning if not found.
    DataSet* CheckForSet( std::string const& ) const;
    /// Get DataSet matching specified attibutes.
    DataSet* CheckForSet( std::string const&, int, std::string const&) const;
    /// Get multiple DataSets matching specified argument.
    DataSetList GetMultipleSets( std::string const& ) const;
    /// Select multiple sets, no warning if none found.
    DataSetList SelectSets( std::string const& ) const;
    /// Generate name based on given default and # of DataSets.
    std::string GenerateDefaultName(std::string const&) const;
    /// Add or append to string DataSet
    DataSet* AddOrAppendSet(std::string const&, int, std::string const&,
                            std::vector<std::string> const&);
    /// Add or append to DataSet
    DataSet* AddOrAppendSet(std::string const&, int, std::string const&,
                            std::vector<double> const&, std::vector<double> const&);
    /// Add DataSet to list with name, or default name if not specified.
    DataSet* AddSet( DataSet::DataType, std::string const&, const char*);
    /// Add DataSet to list with name and index.
    DataSet* AddSetIdx( DataSet::DataType, std::string const&, int);
    /// Add DataSet to list with name and aspect.
    DataSet* AddSetAspect( DataSet::DataType, std::string const&, std::string const&);
    /// Add DataSet to list with name, idx, and aspect.
    DataSet* AddSetIdxAspect( DataSet::DataType, std::string const&, int, std::string const&);
    /// Add DataSet to list with name, idx, aspect, and legend.
    DataSet* AddSetIdxAspect( DataSet::DataType, std::string const&, int, std::string const&,
                              std::string const&);
    /// Add already set up DataSet to list.
    int AddSet( DataSet* );
    /// Add a copy of the DataSet to the list; memory for DataSet will not be freed.
    void AddCopyOfSet(DataSet*);
    /// Print info on DataSets in the list
    void List() const;
#   ifdef MPI
    /// Call sync for DataSets in the list (MPI only)
    void SynchronizeData();
#   endif
    /// Find next set of specified type with given name.
    DataSet* FindSetOfType(std::string const&, DataSet::DataType) const;
    /// Find COORDS DataSet or create default COORDS DataSet.
    DataSet* FindCoordsSet(std::string const&);
    /// Get reference frame DataSet from name/tag
    DataSet* GetReferenceFrame(std::string const&) const;
    /// reference arg help text
    static const char* RefArgs;
    /// Get reference frame DataSet from args
    ReferenceFrame GetReferenceFrame(ArgList&) const;
    /// List all reference frames.
    void ListReferenceFrames() const;
  private:
    /// Separate input string into DataSet args.
    static std::string ParseArgString(std::string const&, std::string&, std::string&);
    /// Warn if DataSet not found but may be pending.
    inline void PendingWarning() const;

    typedef std::vector<DataSet*> DataListType;
    /// DataSet debug level
    int debug_;
    /// Ensemble member number
    int ensembleNum_;
    /// True if list contains copies that should not be freed in destructor.
    bool hasCopies_;
    /// True if Actions will generate DataSets in the future.
    bool dataSetsPending_;
    /// List of DataSets
    DataListType DataList_;
    /// Hold descriptions and allocators for all DataSet types.
    struct DataToken {
      const char* Description;
      DataSet::AllocatorType Alloc;
    };
    static const DataToken DataArray[];
    typedef const DataToken* TokenPtr;
};
#endif
