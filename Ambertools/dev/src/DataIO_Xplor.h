#ifndef INC_DATAIO_XPLOR_H
#define INC_DATAIO_XPLOR_H
#include "DataIO.h"
/// Write Xplor format data files.
class DataIO_Xplor : public DataIO {
  public:
    DataIO_Xplor() : DataIO(false,false,true) {} // Valid for 3D only
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Xplor(); }
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(std::string const&, DataSetList&, std::string const&);
    int processWriteArgs(ArgList&)                 { return 0; }
    int WriteData(std::string const&,DataSetList const&)         { return 1; }
    int WriteData2D(std::string const&, DataSetList const&)      { return 1; }
    int WriteData3D(std::string const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    int WriteSet3D(DataSet const&, CpptrajFile&);
    std::string title_;
    std::string remark_;
};
#endif
