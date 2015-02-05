#include "DataSet_MatrixDbl.h"
void DataSet_MatrixDbl::Write2D(CpptrajFile& outfile, int xIn, int yIn) const {
  size_t x = (size_t)xIn;
  size_t y = (size_t)yIn;
  if ( xIn < 0 || yIn < 0 || x >= mat_.Ncols() || y >= mat_.Nrows() )
    outfile.Printf(data_format_, 0.0);
  else 
    outfile.Printf(data_format_, mat_.element(x,y));
}

double* DataSet_MatrixDbl::MatrixArray() const {
  double* matOut = new double[ mat_.size() ];
  std::copy( mat_.Ptr(), mat_.Ptr() + mat_.size(), matOut );
  return matOut;
}
