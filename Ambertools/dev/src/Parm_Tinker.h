#ifndef INC_PARM_TINKER_H
#define INC_PARM_TINKER_H
#include "ParmIO.h"
class Parm_Tinker : public ParmIO {
  public :
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_Tinker(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&) { return 1; }
    void SetDebug(int) {}
    int processWriteArgs(ArgList&) { return 0; }
};
#endif
