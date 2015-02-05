#include "DataSet_GridFlt.h"
void DataSet_GridFlt::Write3D(CpptrajFile& outfile, int xIn, int yIn, int zIn) const {
  size_t x = xIn;
  size_t y = yIn;
  size_t z = zIn;
  if ( xIn < 0 || yIn < 0 || zIn < 0 ||
       x >= grid_.NX() || y >= grid_.NY() || z >= grid_.NZ() )
    outfile.Printf(data_format_, 0.0);
  else
    outfile.Printf(data_format_, grid_.element(x,y,z));
}
