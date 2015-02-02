#ifndef INC_PARM_CHARMMPSF_H
#define INC_PARM_CHARMMPSF_H
#include "ParmIO.h"
class Parm_CharmmPsf : public ParmIO {
  public :
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_CharmmPsf(); }
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(std::string const&, Topology&);
    int WriteParm(std::string const&, Topology const&);
    void SetDebug(int) {}
    int processWriteArgs(ArgList&) { return 0; }
};
#endif
