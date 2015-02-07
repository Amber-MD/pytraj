#include "DataIO.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_Mesh.h"

// DataIO::CheckValidFor()
bool DataIO::CheckValidFor( DataSet const& dataIn ) const {
  if (valid1d_ && dataIn.Ndim() == 1) return true; 
  if (valid2d_ && dataIn.Ndim() == 2) return true; 
  if (valid3d_ && dataIn.Ndim() == 3) return true;
  for (std::vector<DataSet::DataType>::const_iterator dt = valid_.begin(); 
                                                      dt != valid_.end(); ++dt)
    if (dataIn.Type() == *dt) return true;
  return false;
}

// DataIO::SetupCoordFormat()
std::string DataIO::SetupCoordFormat(size_t maxFrames, Dimension const& dim, 
                                     int default_width, int default_precision)
{
  int col_precision = default_precision;
  // Determine maximum coordinate.
  double maxCoord = (dim.Step() * (double)maxFrames) + dim.Min();
  // Determine character width necessary to hold largest coordinate.
  int col_width = DigitWidth( (long int)maxCoord );
  // Check if the precision is enough to support the step size.
  if (dim.Step() < 1.0) {
    int prec_exp_width = FloatWidth( dim.Step() );
    if (prec_exp_width > col_precision)
      col_precision = prec_exp_width;
  }
  // If the width for the column plus the characters needed for precision
  // (plus 1 for decimal point) would be greated than default_width, increment 
  // the column width by (precision+1).
  if (col_precision != 0) {
    int precision_width = col_width + col_precision + 1;
    if (precision_width > default_width) col_width = precision_width;
  }
  // Default width for column is at least default_width.
  if (col_width < default_width) col_width = default_width;
  // Set column data format string, left-aligned (no leading space).
  return SetDoubleFormatString( col_width, col_precision, 0 );
}
