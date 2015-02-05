#ifndef INC_DATAIO_VECTRAJ_H
#define INC_DATAIO_VECTRAJ_H
#include "DataIO.h"
#include "TrajectoryFile.h"
/// Write vector data to pseudo trajectory.
class DataIO_VecTraj : public DataIO {
  public:
    DataIO_VecTraj();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new DataIO_VecTraj(); }
    static void WriteHelp();
    int processReadArgs(ArgList&) { return 0; }
    int ReadData(std::string const&,DataSetList&,std::string const&) {return 1;}
    int processWriteArgs(ArgList&);
    int WriteData(std::string const&,DataSetList const&);
    int WriteData2D(std::string const&, DataSetList const&) { return 1; }
    int WriteData3D(std::string const&, DataSetList const&) { return 1; }
    bool ID_DataFormat(CpptrajFile&) { return false; }
  private:
    TrajectoryFile::TrajFormatType trajoutFmt_;
    std::string parmoutName_;
};
#endif
