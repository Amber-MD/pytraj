#ifndef INC_ACTION_UNWRAP_H
#define INC_ACTION_UNWRAP_H
#include "Action.h"
#include "ImageTypes.h"
class Action_Unwrap : public Action {
  public:
    Action_Unwrap();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Unwrap(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    Image::PairType imageList_;
    Image::Mode imageMode_;
    AtomMask mask_;
    Frame RefFrame_;
    Topology* RefParm_;
    bool orthogonal_;
    bool center_;
};
#endif
