#ifndef INC_DATAIO_H
#define INC_DATAIO_H
#include "ArgList.h"
#include "DataSetList.h"
#include "CpptrajFile.h"
#include "BaseIOtype.h"
#include "DataSet_double.h"
/// Base class that all DataIO objects inherit from.
/** Reading in data occurs as soon as ReadData is called, so any read-specific
  * arguments are passed in with ReadData. Writing data can occur any time
  * after a DataFile is set up for write, so write arguments are processed
  * separately with processWriteArgs.
  */
class DataIO : public BaseIOtype {
  public:
    DataIO() : debug_(0), valid1d_(false), valid2d_(false), valid3d_(false) {}
    DataIO(bool v1, bool v2, bool v3) :
               valid1d_(v1), valid2d_(v2), valid3d_(v3) {}
    virtual ~DataIO() {}
    // ----- Inherited Functions -----------------
    virtual int processReadArgs(ArgList&) = 0;
    virtual int ReadData(std::string const&,DataSetList&,std::string const&) = 0;
    virtual int processWriteArgs(ArgList&) = 0;
    virtual int WriteData(std::string const&, DataSetList const&) = 0;
    virtual int WriteData2D(std::string const&, DataSetList const&) = 0;
    virtual int WriteData3D(std::string const&, DataSetList const&) = 0;
    virtual bool ID_DataFormat(CpptrajFile&) = 0; // TODO: -> BaseIOtype?
    /// \return True if this DataIO valid for given DataSet
    bool CheckValidFor(DataSet const&) const;
    void SetDebug(int d) { debug_ = d; }
  protected:
    // TODO: Move this to DataSet?
    static std::string SetupCoordFormat(size_t, Dimension const&, int, int);
    /// Indicate this DataIO is valid for given DataSet type
    void SetValid(DataSet::DataType t) { valid_.push_back( t ); }
    int debug_;
  private:
    std::vector<DataSet::DataType> valid_; ///< Data sets for which DataIO is valid writer.
    bool valid1d_; ///< Valid for all 1D data sets.
    bool valid2d_; ///< Valid for all 2D data sets.
    bool valid3d_; ///< Valid for all 3D data sets.
};
#endif
