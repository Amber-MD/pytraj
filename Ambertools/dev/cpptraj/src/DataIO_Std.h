#ifndef INC_DATAIO_STD_H
#define INC_DATAIO_STD_H
#include "DataIO.h"
/// Read/write standard data files.
class DataIO_Std : public DataIO {
  public:
    DataIO_Std();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_Std(); }
    static void ReadHelp();
    static void WriteHelp();
    int processReadArgs(ArgList&);
    int ReadData(std::string const&,DataSetList&,std::string const&);
    int processWriteArgs(ArgList&);
    int WriteData(std::string const&,DataSetList const&);
    int WriteData2D(std::string const&, DataSetList const&);
    int WriteData3D(std::string const&, DataSetList const&);
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    static const char* SEPARATORS;
    int Read_1D(std::string const&,DataSetList&,std::string const&);
    int Read_2D(std::string const&,DataSetList&,std::string const&);
    int Read_3D(std::string const&,DataSetList&,std::string const&);
    int Read_Vector(std::string const&,DataSetList&,std::string const&);
    static void WriteNameToBuffer(CpptrajFile&, std::string const&, int,  bool);
    int WriteDataNormal(CpptrajFile&,DataSetList const&);
    int WriteDataInverted(CpptrajFile&,DataSetList const&);
    int WriteSet2D(DataSet const&, CpptrajFile&);
    int WriteSet3D(DataSet const&, CpptrajFile&);
    enum modeType {READ1D=0, READ2D, READVEC};
    modeType mode_;
    int indexcol_;
    bool isInverted_;  ///< For 1D writes invert X/Y.
    bool hasXcolumn_;
    bool writeHeader_;
    bool square2d_;
};
#endif
