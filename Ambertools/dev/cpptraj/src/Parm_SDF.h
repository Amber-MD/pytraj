#ifndef INC_PARM_SDF_H
#define INC_PARM_SDF_H
#include "ParmIO.h"
class Parm_SDF : public ParmIO {
  public :
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_SDF(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&) { return 1; }
    void SetDebug(int) {}
    int processWriteArgs(ArgList&) { return 0; }
};
#endif
