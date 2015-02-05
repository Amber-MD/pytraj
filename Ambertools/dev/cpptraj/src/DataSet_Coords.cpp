#include "DataSet_Coords.h"

Frame DataSet_Coords::AllocateFrame() const {
  Frame f;
  f.SetupFrameV( top_.Atoms(), hasVel_, 0 );
  return f;
}

void DataSet_Coords::SetTopology(Topology const& topIn) {
  top_ = topIn;
  numCrd_ = top_.Natom() * 3;
  if (top_.ParmBox().HasBox())
    numBoxCrd_ = 6;
  else
    numBoxCrd_ = 0;
  hasVel_ = top_.HasVelInfo();
}
