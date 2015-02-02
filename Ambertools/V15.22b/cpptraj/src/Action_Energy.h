#ifndef INC_ACTION_ENERGY_H
#define INC_ACTION_ENERGY_H
#include "Action.h"
#include "Energy.h"
/// Calculate energy 
class Action_Energy: public Action {
  public:
    Action_Energy();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Energy(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print();
    /// Corresponds to data sets.
    enum Etype { BOND = 0, ANGLE, DIHEDRAL, V14, Q14, VDW, ELEC, TOTAL};
    /// Add energy data set of specified type.
    int AddSet(Etype, DataSetList*, DataFile*, std::string const&);
    /// Corresponds to calculations.
    enum CalcType { BND, ANG, DIH, N14, NBD };
    std::vector<DataSet*> Energy_; ///< Hold output data sets
    std::vector<CalcType> Ecalcs_; ///< Hold which calcs to perform
    typedef std::vector<CalcType>::const_iterator calc_it;
    Topology* currentParm_;        ///< Hold current topology
    AtomMask Mask1_;               ///< Char mask for all but NB calc
    AtomMask Imask_;               ///< Int mask for NB calc
    Energy_Amber ENE_;             ///< Energy calc class.
};
#endif
