// -*- mode: c++; -*-

#ifndef INC_ACTION_PAIRDIST_H
#define INC_ACTION_PAIRDIST_H

#include "Action.h"
#include "ImagedAction.h"
#include "OnlineVarT.h"



/** \author Hannes H. Loeffler
  */

class Action_PairDist : public Action, ImagedAction {
 public:
  Action_PairDist();

  static DispatchObject* Alloc() {
    return (DispatchObject*)new Action_PairDist();
  }

  static void Help();

 private:
  Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
		       DataFileList*, int);
  Action::RetType Setup(Topology*, Topology**);
  Action::RetType DoAction(int, Frame*, Frame**);
  void Print();

  CpptrajFile output_;

  AtomMask mask1_;
  AtomMask mask2_;

  double delta_;		// resolution

  std::vector<Stats<double> > histogram_;
  unsigned long maxbin_;

  bool same_mask_;
  unsigned long ub1_;
  unsigned long ub2_;
};
#endif
