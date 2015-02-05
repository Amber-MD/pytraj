#ifndef INC_ACTION_PUCKER_H
#define INC_ACTION_PUCKER_H
// Class: Action_Pucker
/// Calculate the ring pucker given 5 atom masks.
#include "Action.h"
class Action_Pucker: public Action {
  public:
    Action_Pucker();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Pucker(); }
    static void Help();
  private:
    DataSet *pucker_;
    DataSet* amplitude_;
    DataSet* theta_;
    std::vector<AtomMask> Masks_;
    std::vector<Vec3> AX_;
    enum PmethodType { ALTONA=0, CREMER };
    PmethodType puckerMethod_;
    bool useMass_;
    bool range360_;
    double offset_;

    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
};
#endif  
