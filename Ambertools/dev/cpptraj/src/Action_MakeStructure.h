#ifndef INC_ACTION_MAKESTRUCTURE_H
#define INC_ACTION_MAKESTRUCTURE_H
#include "Action.h"
#include "DihedralSearch.h"
class Action_MakeStructure : public Action {
  public:
    Action_MakeStructure();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_MakeStructure(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, DataSetList*, DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}
    Topology* CurrentParm_; // DEBUG
    int debug_;
    /// Hold secondary structure/turn/single dihedral types.
    class SS_TYPE {
      public: 
        SS_TYPE() {}
        SS_TYPE(double ph,double ps,double ph2,double ps2,int t,std::string const& n) :
          phi(ph), psi(ps), phi2(ph2), psi2(ps2), isTurn(t), type_arg(n) {}
        bool empty() { return (isTurn == -1); }
        double phi, psi, phi2, psi2; // Angle(s) FIXME: Should this be an array?
        int isTurn; // 0=phi/psi, 1=Turn (phi/psi/phi2/psi2), 2=Single Dihedral(phi)
        std::string type_arg;
    };
    std::vector<SS_TYPE> SS;
    /// Determine if SS type has already been defined.
    int FindSStype(std::string const&);
    /// Hold what angles will be applied to which residues.
    class SecStructHolder {
      public:
        SecStructHolder() : sstype_idx(-1) {}
        SecStructHolder(std::string const& rangearg, int typeidx) :
          resRange(rangearg, -1), sstype_idx(typeidx) {}
        Range resRange;                ///< Residues to set phi/psi for.
        DihedralSearch dihSearch_;     ///< Used to search for specified dihedrals
        std::vector<AtomMask> Rmasks_; ///< Masks of atoms to hold fixed during rotation.
        std::vector<float> thetas_;    ///< theta for each dihedral.
        int sstype_idx;                ///< Pointer to corresponding SS_TYPE.
    };
    std::vector<SecStructHolder> secstruct_;
};
#endif
