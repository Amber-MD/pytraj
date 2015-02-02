// -*- mode: c++; -*-

#ifndef INC_ACTION_DENSITY_H
#define INC_ACTION_DENSITY_H

#include <vector>
#include <map>

#include "Action.h"
#include "ImagedAction.h"
#include "OnlineVarT.h"



/** \author Hannes H. Loeffler
  */

#define ROUTINE_VERSION_STRING "1.0.1"

typedef StatsMap<long,double> statmap;

class Action_Density : public Action, ImagedAction {
public:
  Action_Density();

  static DispatchObject* Alloc() {
    return (DispatchObject*) new Action_Density();
  }
    
  static void Help();

private:
  Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
		       DataFileList*, int);
  Action::RetType Setup(Topology*, Topology**);
  Action::RetType DoAction(int, Frame*, Frame**);
  void Print();
  void Output(long, long, std::vector<statmap>&);

  static const std::string emptystring;

  enum DirectionType {DX = 0, DY, DZ};
  enum PropertyType {NUMBER = 0, MASS, CHARGE, ELECTRON};

  DirectionType axis_;
  DirectionType area_coord_[2];
  PropertyType property_;

  double delta_;
  Stats<double> area_;

  CpptrajFile output_;

  std::vector<AtomMask> masks_;

  // std::unordered_map may be better but it is C++11, some STL's may have
  // hash_map but not same number of params,
  // two separate maps are used to ensure we store negative and positive zeros,
  // i.e. to properly bin ]-1,0[ and [0,1[
  std::vector<statmap> minus_histograms_, plus_histograms_;

  std::vector<std::vector<double> > properties_;
};
#endif    
