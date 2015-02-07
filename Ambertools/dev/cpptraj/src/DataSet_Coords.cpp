#include "DataSet_Coords.h"
#include "CpptrajStdio.h"

Frame DataSet_Coords::AllocateFrame() const {
  Frame f;
  f.SetupFrameV( top_.Atoms(), top_.ParmCoordInfo() );
  return f;
}

void DataSet_Coords::SetTopology(Topology const& topIn) {
  top_ = topIn;
  numCrd_ = top_.Natom() * 3;
  if (top_.ParmBox().HasBox())
    numBoxCrd_ = 6;
  else
    numBoxCrd_ = 0;
  hasVel_ = top_.ParmCoordInfo().HasVel(); // TODO: Remove, use CoordinateInfo
  // FIXME: The COORDS DataSet cannot store things like rep dims,
  //        times, or temps. Remove these from the CoordinateInfo
  //        and warn. This should be in DataSet_Coords_CRD itself.
  if (Type() == DataSet::COORDS) {
    CoordinateInfo const& pInfo = top_.ParmCoordInfo();
    if (pInfo.ReplicaDimensions().Ndims() > 0)
      mprintf("Warning: COORDS data sets do not store replica dimensions.\n");
    if (pInfo.HasTemp())
      mprintf("Warning: COORDS data sets do not store temperatures.\n");
    if (pInfo.HasTime())
      mprintf("Warning: COORDS data sets do not store times.\n");
    top_.SetParmCoordInfo(CoordinateInfo(pInfo.TrajBox(),pInfo.HasVel(),false,false));
  }
}
