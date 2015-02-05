#include <cmath> // regression
#include "DataSet_Mesh.h"
#include "CpptrajStdio.h"
#include "Constants.h" // regression

// CONSTRUCTOR - Create X mesh
DataSet_Mesh::DataSet_Mesh(int sizeIn, double ti, double tf) :
  DataSet_1D(XYMESH, 12, 4)
{
  CalculateMeshX(sizeIn, ti, tf);
}

// DataSet_Mesh::Allocate1D()
int DataSet_Mesh::Allocate1D( size_t sizeIn ) {
  mesh_x_.reserve( sizeIn );
  mesh_y_.reserve( sizeIn );
  return 0;
}

// DataSet_Mesh::WriteBuffer()
void DataSet_Mesh::WriteBuffer(CpptrajFile &cbuffer, size_t frame) const {
  if (frame >= mesh_x_.size())
    cbuffer.Printf(data_format_, 0.0);
  else
    cbuffer.Printf(data_format_, mesh_y_[frame]);
}

// -----------------------------------------------------------------------------
void DataSet_Mesh::Append(std::vector<double> const& xIn,
                          std::vector<double> const& yIn)
{
  if (xIn.empty() || xIn.size() != yIn.size()) return;
  size_t oldsize = Size();
  mesh_x_.resize( oldsize + xIn.size() );
  mesh_y_.resize( oldsize + yIn.size() );
  std::copy( xIn.begin(), xIn.end(), mesh_x_.begin() + oldsize );
  std::copy( yIn.begin(), yIn.end(), mesh_y_.begin() + oldsize );
}

// DataSet_Mesh::CalculateMeshX()
void DataSet_Mesh::CalculateMeshX(int sizeIn, double ti, double tf) {
  mesh_x_.resize( sizeIn, 0 );
  mesh_y_.resize( sizeIn, 0 );
  double s = (ti + tf)/2;
  double d = (tf - ti)/2;
  for (int i = 0; i < sizeIn; i++)
    mesh_x_[i] = s + d*((double) (2*i + 1 - sizeIn)/(sizeIn - 1));
}

// DataSet_Mesh::SetMeshXY()
int DataSet_Mesh::SetMeshXY(DataSet_1D const& dsIn) {
  mesh_x_.resize( dsIn.Size() );
  mesh_y_.resize( dsIn.Size() );
  for (unsigned int i = 0; i < dsIn.Size(); i++) {
    mesh_x_[i] = dsIn.Xcrd(i);
    mesh_y_[i] = dsIn.Dval(i);
  }
  return 0;
}

// ---------- Integration routines ---------------------------------------------
// DataSet_Mesh::Integrate_Trapezoid()
double DataSet_Mesh::Integrate_Trapezoid( DataSet_Mesh& sumOut ) const {
  double sum = 0.0;
  int mesh_size = (int)mesh_x_.size();
  if (mesh_size < 2) return 0.0;
  // Give output data set the same X mesh
  sumOut.mesh_x_ = mesh_x_;
  sumOut.mesh_y_.resize( mesh_x_.size() );
  sumOut.mesh_y_[0] = 0.0;
  for (int i = 1; i < mesh_size; i++) {
      double b_minus_a = (mesh_x_[i] - mesh_x_[i - 1]);
      sum += (b_minus_a * (mesh_y_[i - 1] + mesh_y_[i]) * 0.5);
      sumOut.mesh_y_[i] = sum;
  }
  return sum;
}

// DataSet_Mesh::Integrate_Trapezoid()
double DataSet_Mesh::Integrate_Trapezoid() const {
  double sum = 0.0;
  int mesh_size = (int)mesh_x_.size();
  if (mesh_size < 2) return 0.0;
  for (int i = 1; i < mesh_size; i++) {
      double b_minus_a = (mesh_x_[i] - mesh_x_[i - 1]);
      sum += (b_minus_a * (mesh_y_[i - 1] + mesh_y_[i]) * 0.5);
  }
  return sum;
}

// ---------- Cubic Spline Routines --------------------------------------------
// DataSet_Mesh::cubicSpline_coeff()
/** Given a set of x and y values of size n, compute the b, c, and d
  * coefficients for n interpolating cubic splines of form:
  *
  *   Si(t) = y[i] + b[i]*dt + c[i]*dt^2 + d[i]*dt^3 
  *
  * where dt = t - x[i]. Si(t), and first and second derivatives Si'(t)
  * and Si"(t) must be continuous over interval dt, and both derivatives 
  * for adjacent points must be equal so that adjacent line segments become 
  * continuous.
  *
  *   Si'(t) = b[i] + 2*c[i]*dt + 3*d[i]*dt^2
  *   Si"(t) = 2*c[i] + 6*d[i]*dt
  *
  * \param x Input X values
  * \param y Corresponding Y values
  */
void DataSet_Mesh::cubicSpline_coeff(std::vector<double> const& x, std::vector<double> const& y) 
{
  int n = (int)x.size();

  b.resize(n, 0.0);
  c.resize(n, 0.0);
  d.resize(n, 0.0);

  int n_minus1 = n - 1;

  if ( n > 2 ) {
    // Generate Tri-diagonal matrix
    d[0] = x[1] - x[0];
    c[1] = (y[1] - y[0]) / d[0];
    for (int i = 1; i < n_minus1; i++) {
      d[i] = x[i + 1] - x[i];
      b[i] = 2.0 * (d[i - 1] + d[i]);
      c[i+1] = (y[i + 1] - y[i]) / d[i];
      c[i] = c[i+1] - c[i];
      //mprintf("TRIDBG: %i b=%lf c=%lf d=%lf\n",i,b[i],c[i],d[i]);
    }
    
    // Set up boundary  conditions
    b[0]        = -d[0];
    b[n_minus1] = -d[n - 2];
    c[0]        = 0.0;
    c[n_minus1] = 0.0;
    if (n > 3) {
      c[0]        = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0]);
      c[n_minus1] = c[n - 2] / (x[n_minus1] - x[n - 3]) - c[n - 3] / (x[n - 2] - x[n - 4]);
      c[0]        = c[0] * d[0] * d[0] / (x[3] - x[0]);
      c[n_minus1] = -c[n_minus1] * d[n - 2] * d[n - 2] / (x[n_minus1] - x[n - 4]);
    }

    // Forward elimination
    for (int i = 1; i < n; i++) {
        double t = d[i - 1] / b[i - 1];
        b[i]     = b[i] - t * d[i - 1];
        c[i]     = c[i] - t * c[i - 1];
        //mprintf("FWDDBG: %i b=%lf c=%lf t=%lf\n",i,b[i],c[i],t);
    }

    // Back substitution
    c[n_minus1] = c[n_minus1] / b[n_minus1];
    for (int i = n - 2; i > -1; i--) {
        c[i] = (c[i] - d[i] * c[i + 1]) / b[i];
        //mprintf("BAKDBG: %i c=%lf\n",i,c[i]);
    }

    // Calculate the polynomial coefficients
    b[n_minus1] = (y[n_minus1] - y[n - 2]) / d[n - 2] + d[n - 2] * (c[n - 2] + 2.0 * c[n_minus1]);
    for (int i = 0; i < n_minus1; i++) {
        b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i]);
        d[i] = (c[i + 1] - c[i]) / d[i];
        c[i] = 3.0 * c[i];
        //mprintf("POLYDBG: %i b=%lf c=%lf d=%lf\n",i,b[i],c[i],d[i]);
    }
    c[n_minus1] = 3.0 * c[n_minus1];
    d[n_minus1] = d[n - 2];

  // Special case; n == 2
  } else {
    b[0] = (y[1] - y[0]) / (x[1] - x[0]);
    c[0] = 0.0;
    d[0] = 0.0;
    b[1] = b[0];
    c[1] = 0.0; 
    d[1] = 0.0;
  }
}

// DataSet_Mesh::cubicSpline_eval() 
/** Evaluate cubic spline function with pre-calcd coefficients in b, c, and
  * d from coordinates x/y for all points in mesh.
  * \param x Input X coordinates
  * \param y Input Y coordinates
  */
void DataSet_Mesh::cubicSpline_eval(std::vector<double> const& x, std::vector<double> const& y)
{
  int xidx;  
  int n = (int)x.size();
  int mesh_size = (int)mesh_x_.size();

  for (int uidx = 0; uidx < mesh_size; uidx++) {
    double U = mesh_x_[uidx];
    // Search for U in x
    if (U < x[0])
      xidx = 0;
    else if (U > x[n-1])
      xidx = n - 1;
    else {
      int i0 = 0;
      int i1 = n - 1;
      while (i0 <= i1) {
        xidx = (i0 + i1) / 2;
        if ( U < x[xidx] )
          i1 = xidx - 1;
        else if ( U > x[xidx+1] )
          i0 = xidx + 1;
        else
          break;
      }
    }
    // Evaluate v for this u
    double dx = U - x[xidx];
    mesh_y_[uidx] = y[xidx] + dx*(b[xidx] + dx*(c[xidx] + dx*d[xidx])); 
  }
}

// DataSet_Mesh::SetSplinedMeshY()
/** Assumes mesh X values already set with CalculateMeshX. */
int DataSet_Mesh::SetSplinedMeshY(std::vector<double> const& x, std::vector<double> const& y) {
  if (x.size() != y.size()) {
    mprinterr("Error: X size (%u) != Y size (%u)\n", x.size(), y.size());
    return 1;
  }
  // No point if 1 or less values
  if (x.size() < 2) {
    mprinterr("Error: Requires > 1 values (%u specified).\n", x.size());
    return 1;
  }
  cubicSpline_coeff(x, y);
  cubicSpline_eval(x, y);
  return 0;
}

// DataSet_Mesh::SetSplinedMesh()
/** Assumes mesh X values already set with CalculateMeshX. */
int DataSet_Mesh::SetSplinedMesh(DataSet_1D const& dsIn)
{
  if (dsIn.Size() < 2) {
    mprinterr("Error: Requires > 1 values (%u specified).\n", dsIn.Size());
    return 1;
  }
  // Create X and Y values for dsIn
  std::vector<double> x, y;
  x.reserve( dsIn.Size() );
  y.reserve( dsIn.Size() );
  for (int i = 0; i < (int)dsIn.Size(); i++) {
    x.push_back( dsIn.Xcrd( i ) );
    y.push_back( dsIn.Dval( i ) );
  }
  cubicSpline_coeff(x, y);
  cubicSpline_eval(x, y);
  return 0;
}

// ---------- Linear Regression ------------------------------------------------
/** This code (especially the error analysis) was adapted from grace 5.1.22
  * fit.c:linear_regression().
  */
int DataSet_Mesh::LinearRegression( double& slope, double& intercept, 
                                    double& correl, bool silent ) const
{
  if (mesh_x_.size() < 2) return 1;
  double mesh_size = (double)mesh_x_.size();
  // Averages
  double xavg = 0.0, yavg = 0.0;
  for (unsigned int i = 0; i < mesh_x_.size(); i++) {
    xavg += mesh_x_[i];
    yavg += mesh_y_[i];
  }
  xavg /= mesh_size;
  yavg /= mesh_size;
  // Sums of squares
  double sxx = 0.0, sxy = 0.0, syy = 0.0;
  for (unsigned int i = 0; i < mesh_x_.size(); i++) {
    double xdiff = mesh_x_[i] - xavg;
    double ydiff = mesh_y_[i] - yavg;
    sxx += (xdiff * xdiff);
    sxy += (xdiff * ydiff);
    syy += (ydiff * ydiff);
  }
  // Standard deviation, covariance
  double xsd = sqrt( sxx / (mesh_size - 1.0) );
  double ysd = sqrt( syy / (mesh_size - 1.0) );
  if (xsd < Constants::SMALL || ysd < Constants::SMALL) {
    mprinterr("Error: All values of x or y are the same (SD cannot be zero).\n");
    return 1;
  }
  double covariance = sxy / (mesh_size - 1.0);
         correl = covariance / (xsd * ysd);
         slope = sxy / sxx;
         intercept = yavg - slope * xavg;
  if (!silent) {
    mprintf("\tData points= %u\n"
            "\t<X>= %g\n\t<Y>= %g\n"
            "\tSDx= %g\n\tSDy= %g\n"
            "\tCorrelation coefficient= %g\n"
            "\tSlope= %g\n", mesh_x_.size(),
            xavg, yavg, xsd, ysd, correl, slope);
  }
  // Case N==2, no need for error analysis.
  if (mesh_x_.size() == 2) {
    slope = (mesh_y_[1] - mesh_y_[0]) / (mesh_x_[1] - mesh_x_[0]);
    intercept = mesh_y_[0] - slope * mesh_x_[0];
    if (!silent) mprintf("\tIntercept= %g\n", intercept);
    return 0;
  } 
  // Error analysis
  double residualSumSq = syy - slope * sxy;
  double residualMeanSq = residualSumSq / (mesh_size - 2.0);
  //double stdErrRegression = sqrt( residualMeanSq );
  double stdErrIntercept = sqrt( residualMeanSq * (1.0 / mesh_size + xavg * xavg / sxx) );
  double stdErrSlope = sqrt( residualMeanSq / sxx );
  double sumSqRegression = syy - residualSumSq;
  double Fval = sumSqRegression / residualMeanSq;
  //double R2 = sumSqRegression / syy;
  if (!silent) {
    mprintf("\tStandard error of slope= %g\n"
            "\tt - value for slope= %g\n"
            "\tIntercept= %g\n"
            "\tStandard Error of intercept= %g\n"
            "\tt - value for intercept= %g\n",
            stdErrSlope, slope / stdErrSlope, 
            intercept, stdErrIntercept, intercept / stdErrIntercept);

    mprintf("\tVariance analysis:\n\t%-10s %5s %14s %14s %14s\n",
            "Source", "d.f", "Sum of squares", "Mean square", "F");
    mprintf("\t%-10s %5d %14.7g %14.7g %14.7g\n", "Regression",
            1, sumSqRegression, sumSqRegression, Fval);
    mprintf("\t%-10s %5u %14.7g %14.7g\n", "Residual",
            mesh_x_.size() - 2, residualSumSq, residualMeanSq);
    mprintf("\t%-10s %5u %14.7g\n", "Total",  mesh_x_.size() - 1, syy);
  }
  return 0;
}
