#include "DataSet_Coords_CRD.h"
#include "CpptrajStdio.h"

int DataSet_Coords_CRD::Allocate1D( size_t sizeIn ) {
  coords_.reserve( sizeIn );
  return 0;
}

void DataSet_Coords_CRD::Info() const {
  size_t sze = (((coords_.size() * (size_t)top_.Natom() * 3UL) + numBoxCrd_) * sizeof(float))
               / 1048576UL;
  if (sze == 0)
    mprintf(" (<1 MB)");
  else
    mprintf(" (%zu MB)", sze);
  if (numBoxCrd_ > 0) mprintf(" Box Coords,");
  if (hasVel_)        mprintf(" Velocities,");
  mprintf(" %i atoms", top_.Natom());
}
