#ifndef INC_ACTION_CREATERESERVOIR_H
#define INC_ACTION_CREATERESERVOIR_H
#include "Action.h"
#include "Traj_AmberNetcdf.h"
#include "DataSet_1D.h"
// Class: Action_CreateReservoir
/// Create a RREMD structure reservoir.
class Action_CreateReservoir : public Action {
  public:
    Action_CreateReservoir();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_CreateReservoir(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
#   ifdef BINTRAJ
    Traj_AmberNetcdf reservoir_;
#   endif
    Topology* original_trajparm_;
    DataSet_1D* ene_;
    DataSet_1D* bin_;
    double reservoirT_;
    int iseed_;
    std::string filename_;
    bool trajIsOpen_;
    bool useVelocity_;
    size_t nframes_;
};
#endif
