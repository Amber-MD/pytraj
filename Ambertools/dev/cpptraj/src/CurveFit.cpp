#ifdef DBG_CURVEFIT
# include <cstdio> // DEBUG
# include <cstdarg> // DEBUG
#endif
#include <cmath> //sqrt
#include <cfloat> // DBL_MIN
#include "CurveFit.h"

// CONSTRUCTOR
CurveFit::CurveFit() : fxn_(0), m_(0), n_(0), errorMessage_(0)
{
# ifdef DBG_CURVEFIT
  // DEBUG
  dbgfile_ = fopen("dbgcurvefit.dat", "wb");
# endif
}

// DESTRUCTOR
#ifdef DBG_CURVEFIT
CurveFit::~CurveFit() { fclose(dbgfile_); }
#endif

void CurveFit::DBGPRINT(const char* format, ...) const {
# ifdef DBG_CURVEFIT
  va_list args;
  va_start(args, format);
  vfprintf(dbgfile_, format, args);
  va_end(args);
# endif
}

const char* CurveFit::Message(int info) {
  switch (info) {
    case 0 : return "Problem with input parameters.";
    case 1 : return "Both actual and predicted relative reductions"
                    " in the sum of squares are at most tolerance.";
    case 2 : return "Relative error between two consecutive iterates"
                    " is at most xtol.";
    case 3 : return "Both actual and predicted relative reductions"
                    " in the sum of squares are at most tolerance"
                    " and relative error between two consecutive iterates"
                    " is at most xtol.";
    case 4 : return "The cosine of the angle between residual and any"
                    " column of the Jacobian is at most gtol in"
                    " absolute value.";
    case 5 : return "Number of calls to function has reached or exceeded max.";
    case 6 : return "ftol is too small. No further reduction in"
                    " the sum of squares is possible.";
    case 7 : return "xtol is too small. no further improvement in"
                    " the approximate solution parameter vector is possible.";
    case 8 : return "gtol is too small. Residual is orthogonal to the"
                    " columns of the Jacobian to machine precision.";
    // The following (9 and 10) are for Statistics() only.
    case 9 : return "Cannot calculate statistics; # elements does not match"
                    " or curve fitting has not been performed.";
    case 10: return "Input set Y values contain zero, cannot calculate RMS"
                    " percent error.";
  }
  return 0;
}

/// Machine precision
// NOTE: from http://cpansearch.perl.org/src/RKOBES/Math-Cephes-0.47/libmd/const.c
const double CurveFit::machine_epsilon = pow(2.0, -53);

// For debugging
void CurveFit::PrintMatrix(const char* name, int ncols, int mrows, CurveFit::Darray const& mat)
const {
  CurveFit::Darray::const_iterator mij = mat.begin();
  for (int n = 0; n < ncols; n++) {
    for (int m = 0; m < mrows; m++)
      DBGPRINT("DEBUG: %s(%i,%i)= %12.6g [%zu]\n", name, m+1, n+1, *(mij++), mij-mat.begin());
  }
}

void CurveFit::PrintVector(const char* name, CurveFit::Darray const& Vec) const {
  DBGPRINT("%s={", name);
  for (CurveFit::Darray::const_iterator v = Vec.begin(); v != Vec.end(); ++v)
    DBGPRINT(" %g", *v);
  DBGPRINT(" }\n");
}

void CurveFit::Pvec_to_Params(Darray& PvecIn) {
  //Params_ = PvecIn;
  static const double epsilon = 1.0e-6;

  for (dsize in = 0; in != n_; in++) {
    if (hasBounds_[in]) {
      double minPmax2 = (Lbound_[in] + Ubound_[in]) / 2.0;
      double maxMmin2 = (Ubound_[in] - Lbound_[in]) / 2.0;
      double t = (PvecIn[in] - minPmax2) / maxMmin2;
      if ( t < -(1.0 - epsilon)) {
        t = -(1.0 - epsilon);
        PvecIn[in] = minPmax2 + maxMmin2 * t;
      }
      if ( t > (1.0 - epsilon)) {
        t = (1.0 - epsilon);
        PvecIn[in] = minPmax2 + maxMmin2 * t;
      }
      Params_[in] = t / (1.0 - fabs(t));
    } else
      Params_[in] = PvecIn[in];
  }
}

void CurveFit::Params_to_Pvec(Darray& PvecIn, Darray const& ParamsIn) const {
  //PvecIn = ParamsIn;
  for (dsize in = 0; in != n_; in++) {
    if (hasBounds_[in]) {
      // map (-inf, inf) to (-1, 1) 
      double t = ParamsIn[in] / (fabs(ParamsIn[in] + 1.0));
      // map (-1,1) to (lower, upper)
      PvecIn[in] = ((Lbound_[in] + Ubound_[in]) / 2.0) +
                   ((Ubound_[in] - Lbound_[in]) / 2.0) * t;
    } else
      PvecIn[in] = ParamsIn[in];
  }
}

/** Evaulate function at given X values and parameters. Calculate residual
  * from given Y values.
  */
void CurveFit::EvaluateFxn(Darray const& Xvals_, Darray const& Yvals_, 
                           Darray const& ParamsIn, Darray& residual)
{
  Params_to_Pvec(fParms_, ParamsIn);
  PrintVector("Param", fParms_);
  fxn_(Xvals_, fParms_, finalY_);
  for (dsize im = 0; im != m_; im++)
  {
    // Residual
    residual[im] = finalY_[im] - Yvals_[im];
  }
  // Apply weights
  for (dsize im = 0; im != Weights_.size(); im++)
    residual[im] *= Weights_[im];
  PrintVector("Residual", residual);
}
  
/** Calculate Jacobian of function from forward-difference approximation.
  * NOTE: The Jacobian is stored as its transpose, i.e. columns are in 
  * rows, to make calculating things like column norms easier. Adapted
  * from subroutine fdjac2_ in lmdif.c from Grace 5.1.22.
  */
// NOTE: ParamsIn is not const so it can be changed when calc. derivative.
void CurveFit::CalcJacobian_ForwardDiff(Darray const& Xvals_, Darray const& Yvals_,
                                        Darray& ParamsIn,
                                        Darray const& residual, Darray& newResidual)
{
  // NOTE: The zero could eventually be a passed-in constant.
  double eps = sqrt( std::max(0.0, machine_epsilon) );
  for (dsize in = 0; in != n_; in++) {
    double param = ParamsIn[in];
    double delta = eps * fabs( param );
    if (delta == 0.0)
      delta = eps;
    ParamsIn[in] = param + delta;
    
    EvaluateFxn(Xvals_, Yvals_, ParamsIn, newResidual);

    ParamsIn[in] = param;

    for (dsize im = 0; im != m_; im++)
      jacobian_[im + in * m_] = (newResidual[im] - residual[im]) / delta;
  }
}
    
/** Calculate normal of given vector within matrix. This routine accumulates
  * the sum of squares in three different sums: small, mid, and large. The
  * sums of the squares for the small and large components are scaled so that
  * no overflows occur. Non-destructive underflows are permitted. The
  * definitions of the small, mid, and large components depend on two
  * constants, rdward and rgiant. The main restrictions on these constants
  * are that rdwarf^2 not underflow and rgiant^2 not overflow.
  * This routine was adapted from the enorm_ subroutine of lmdif.c from
  * Grace 5.1.22, which in turn was adapted from minpack (Argonne National
  * Laboratory. minpack project. March 1980. Burton S. Garbow,
  * Kenneth E. Hillstrom, Jorge J. More).
  */
double CurveFit::VecNorm( Darray::const_iterator const& vBeg, dsize nElt ) {
  // CONSTANTS
  static const double rdwarf = 3.834e-20;
  static const double rgiant = 1.304e19;
  
  double agiant = rgiant / (double)nElt;

  double sum_large = 0.0;
  double sum_mid = 0.0;
  double sum_small = 0.0;
  double lMax = 0.0;
  double sMax = 0.0;
 
  Darray::const_iterator vEnd = vBeg + nElt;
  for (Darray::const_iterator v = vBeg; v != vEnd; ++v)
  {
    double xabs = fabs( *v );
    if (xabs > rdwarf && xabs < agiant)
    {
      // Sum for intermediate components.
      sum_mid += xabs * xabs;
    } else if (xabs <= rdwarf) {
      // Sum for small components
      if (xabs > sMax) {
        double d1 = sMax / xabs;
        sum_small = 1.0 + sum_small * (d1 * d1);
        sMax = xabs;
      } else {
        if (xabs != 0.0) {
          double d1 = xabs / sMax;
          sum_small += d1 * d1;
        }
      }
    } else {
      // Sum for large components
      if (xabs > lMax) {
        double d1 = lMax / xabs;
        sum_large = 1.0 + sum_large * (d1 * d1);
        lMax = xabs;
      } else { // L10
        double d1 = xabs / lMax;
        sum_large += d1 * d1;
      } // L20
    }
  }

  // Calculation of norm
  double ret_val;
  if (sum_large != 0.0) {
    ret_val = lMax * sqrt( sum_large + sum_mid / lMax / lMax );
  } else if (sum_mid != 0.0) {
    if (sum_mid < sMax)
      ret_val = sqrt( sMax * (sum_mid / sMax + sMax * sum_small) );
    else
      ret_val = sqrt( sum_mid * (1.0 + sMax / sum_mid * (sMax * sum_small)));
  } else {
    ret_val = sMax * sqrt(sum_small);
  }
  return ret_val;
}

// CurveFit::ParametersHaveProblems()
bool CurveFit::ParametersHaveProblems(Darray const& Xvals_, Darray const& Yvals_,
                                      Darray const& ParamsIn)
{
  if (ParamsIn.empty() || Xvals_.empty() || Yvals_.empty()) {
    errorMessage_ = "Parameters or coordinates are empty.";
    return true;
  }
  if (Xvals_.size() != Yvals_.size()) {
    errorMessage_ = "Number of X values != number of Y values.";
    return true;
  }
  if ( Xvals_.size() < ParamsIn.size() ) {
    errorMessage_ = "Number of parameters cannot be greater than number of XY values.";
    return true;
  }
  if (hasBounds_.empty())
    hasBounds_.assign(ParamsIn.size(), false);
  else {
    if (hasBounds_.size() != ParamsIn.size() ||
        Ubound_.size() != ParamsIn.size() ||
        Lbound_.size() != ParamsIn.size())
    {
      errorMessage_ = "Number of bounds does not match number of parameters.";
      return true;
    }
    for (dsize in = 0; in != ParamsIn.size(); in++) {
      if (hasBounds_[in]) {
        if ( Lbound_[in] >= Ubound_[in] ) {
          errorMessage_ = "Lower bound must be less than upper bound.";
          return true;
        }
        if (ParamsIn[in] <= Lbound_[in] ||
            ParamsIn[in] >= Ubound_[in])
        {
          errorMessage_ = "Initial parameter not within bounds.";
          return true;
        }
      }
    }
  }
  if (!Weights_.empty() && Weights_.size() != Xvals_.size()) {
    errorMessage_ = "Number of weights does not match number of XY values.";
    return true;
  }
  errorMessage_ = 0;
  return false;
}

/** Perform curve fitting via Levenberg-Marquardt with no bounds/weights. */
int CurveFit::LevenbergMarquardt(FitFunctionType fxnIn, Darray const& Xvals_,
                                 Darray const& Yvals_, Darray& ParamVec,
                                 double tolerance, int maxIterations)
{
  return LevenbergMarquardt(fxnIn, Xvals_, Yvals_, ParamVec, std::vector<bool>(),
                            Darray(), Darray(), Darray(), tolerance, maxIterations);
}

/** Perform curve fitting via Levenberg-Marquard with weights. */
int CurveFit::LevenbergMarquardt(FitFunctionType fxnIn, Darray const& Xvals,
                                 Darray const& Yvals, Darray& ParamVec,
                                 Darray const& weightsIn,
                                 double tolerance, int maxIterations)
{
  return LevenbergMarquardt(fxnIn, Xvals, Yvals, ParamVec, std::vector<bool>(),
                            Darray(), Darray(), weightsIn, tolerance, maxIterations);
}

// -----------------------------------------------------------------------------
/** Perform curve fitting via Levenberg-Marquardt method with optional bounds.
  * \param fxnIn Function to fit.
  * \param Xvals_ target X values (ordinates, M)
  * \param Yvals_ target Y values (coordinates, M)
  * \param ParamVec input parameters(N); at finish contains best esimate of fit parameters
  * \param boundsIn true if parameter has bounds (N)
  * \param lboundIn contain lower bounds for parameter (N)
  * \param uboundIn contain upper bounds for parameter (N)
  * \param weightsIn contain weights for Y values (M)
  * \param tolerance Fit tolerance
  * \param maxIterations Number of iterations to try.
  */
int CurveFit::LevenbergMarquardt(FitFunctionType fxnIn, Darray const& Xvals_,
                                 Darray const& Yvals_, Darray& ParamVec,
                                 std::vector<bool> const& boundsIn,
                                 Darray const& lboundIn, Darray const& uboundIn,
                                 Darray const& weightsIn,
                                 double tolerance, int maxIterations)
{
  int info = 0;
  hasBounds_ = boundsIn;
  Ubound_ = uboundIn;
  Lbound_ = lboundIn;
  Weights_ = weightsIn;
  if (ParametersHaveProblems(Xvals_, Yvals_, ParamVec)) {
    DBGPRINT("Error: %s\n", errorMessage_);
    return info;
  }
  // Set initial parameters
  finalY_ = Yvals_;
  Params_= ParamVec;  // Internal parameter vector
  fParms_ = ParamVec; // Parameters for function evaluation
  Pvec_to_Params( ParamVec );

  fxn_ = fxnIn;
  m_ = Xvals_.size();        // Number of values (rows)
  n_ = Params_.size();       // Number of parameters (cols)
  Darray residual_( m_ );    // Contain Y vals at current Params minus original Y
  // Jacobian matrix/R: m rows by n cols, transposed.
  jacobian_.assign( m_ * n_, 0.0 );

  // For holding diagonal elements for scaling // TODO: Rename
  Darray diag( n_, 0.0 );
  // Workspace arrays
  Darray work1( n_, 0.0 );
  Darray work2( n_, 0.0 );
  Darray newResidual( m_, 0.0 );
  // delta specifies an upper bound on the euclidean norm of d*x
  double delta = 0.0;
  // factor is positive input variable used in determining the initial step
  // bound. This bound is set to the product of factor and the euclidean norm
  // of diag*x if nonzero, or else to factor itself. In most cases factor should
  // lie in the interval (.1,100.). 100. is a generally recommended value.
  // TODO: Make input parameter
  const double factor = 100.0;
  // gtol is a nonnegative input variable. Termination occurs when the cosine 
  // of the angle between residual and any column of the Jacobian is at most 
  // gtol in absolute value. Therefore, gtol measures the orthogonality
  // desired between the function vector and the columns of the Jacobian.
  // TODO: Make input parameter
  const double gtol = 0.0;
  // Smallest possible magnitude
  // TODO: Make Constant?
  const double dwarf = DBL_MIN;
  // Levenberg-Marquard parameter
  double LM_par = 0.0;
  //  xtol is a nonnegative input variable. Termination occurs when the 
  //  relative error between two consecutive iterates is at most xtol.
  //  Therefore, xtol measures the relative error desired in the 
  //  approximate solution.
  // TODO: Make input paramter
  double xtol = 0.0;
  // maxfev is the maxiumum number of allowed function evaluations.
  // NOTE: Included here for complete compatibility with eariler code,
  //       though may not be strictly necessary.
  dsize maxfev = (n_ + 1) * (dsize)maxIterations;

  // Evaluate initial function, obtain residual
  EvaluateFxn( Xvals_, Yvals_, Params_, residual_ );

  // Calculate norm of the residual
  double rnorm = VecNorm( residual_ );
  DBGPRINT("Rnorm= %g\n", rnorm);

  dsize nfev = 1; // Will be set to calls to fxn_. 1 already.

  // MAIN LOOP
  int currentIt = 0;
  while ( currentIt < maxIterations ) {
    DBGPRINT("DEBUG: ----- Iteration %i ------------------------------\n", currentIt+1);
    // Calculate the Jacobian using the forward-difference approximation.
    CalcJacobian_ForwardDiff( Xvals_, Yvals_, Params_, residual_, newResidual );
    PrintMatrix( "Jacobian", n_, m_, jacobian_ );
    nfev += n_;

    // -------------------------------------------
    // | BEGIN QRFAC
    // -------------------------------------------
    // Perform QR factorization of Jacobian via Householder transformations
    // with column pivoting:
    //   J * P = Q * R
    // where J is the Jacobian, P is the permutation matrix, Q is an orthogonal
    // matrix, and R is an upper trapezoidal matrix with diagonal elements of
    // nonincreasing magnitude. The Householder transformation for column k,
    // k = 1,2,...,min(m,n), is of the form:
    //   I = (1/u[k])*u*ut
    // where u has zeros in the first k-1 positions. Adapted from routine qrfac_
    // in Grace 5.1.22 lmdif.c, which in turn is derived from an earlier linpack
    // subroutine.
    DBGPRINT("\nQRFAC ITERATION %i\n", currentIt+1);
    // Rdiag will hold the diagonal elements of R.
    Darray Rdiag(n_, 0.0);
    // Jpvt defines the permutation matrix p such that a*p = q*r. Column j of p
    // is column Jpvt[j] of the identity matrix.
    Iarray Jpvt(n_, 0.0);
    // JcolNorm Will hold the norms of the columns of input Jacobian.
    Darray JcolNorm_(n_, 0.0);
    // Compute the initial column norms and initialize arrays.
    for (dsize in = 0; in != n_; in++) {
      JcolNorm_[in] = VecNorm( jacobian_.begin() + (in * m_), m_ ); 
      Rdiag[in] = JcolNorm_[in];
      work1[in] = Rdiag[in];
      Jpvt[in] = in;
      DBGPRINT("%lu: Rdiag= %12.6g    wa= %12.6g    ipvt= %i\n",
             in+1, Rdiag[in], work1[in], Jpvt[in]+1);
    }

    // Reduce Jacobian to R with Householder transformations.
    dsize min_m_n = std::min( m_, n_ );
    for (dsize in = 0; in != min_m_n; in++) {
      // Bring the column of largest norm into the pivot position.
      dsize kmax = in;
      for (dsize k = in; k != n_; k++) {
        DBGPRINT("\tRdiag[%lu]= %g    Rdiag[%lu]= %g\n", k+1, Rdiag[k], kmax+1, Rdiag[kmax]);
        if (Rdiag[k] > Rdiag[kmax])
          kmax = k;
      }
      DBGPRINT("Elt= %lu    Kmax= %lu\n", in+1, kmax+1);
      if (kmax != in) {
        for (dsize i = 0; i != m_; i++) {
          double temp = jacobian_[i + in * m_];
          jacobian_[i + in   * m_] = jacobian_[i + kmax * m_];
          jacobian_[i + kmax * m_] = temp;
          DBGPRINT("DBG: Swap jac[%lu,%lu] with jac[%lu,%lu]\n", i+1, in+1, i+1, kmax+1);
        }
        Rdiag[kmax] = Rdiag[in];
        work1[kmax] = work1[in];
        dsize k = Jpvt[in];
        Jpvt[in] = Jpvt[kmax];
        Jpvt[kmax] = k;
      }

      // Compute the Householder transformation to reduce the j-th column
      // of Jacobian to a multiple of the j-th unit vector.
      double ajnorm = VecNorm( jacobian_.begin() + (in + in * m_), m_ - in );
      DBGPRINT("#Elt= %lu    mat[%lu]= %g    ajnorm= %g\n",
             m_ - in, in + in * m_, jacobian_[in + in * m_], ajnorm);
      if (ajnorm != 0.0) {
        if (jacobian_[in + in * m_] < 0.0)
          ajnorm = -ajnorm;
        for (dsize im = in; im < m_; im++)
          jacobian_[im + in * m_] /= ajnorm;
        jacobian_[in + in * m_] += 1.0;

        // Apply the transfomation to the remaining columns and update the norms
        dsize in1 = in + 1;
        if (in1 < n_) {
          for (dsize k = in1; k < n_; k++) {
            double sum = 0.0;
            for (dsize i = in; i < m_; i++)
              sum += jacobian_[i + in * m_] * jacobian_[i + k * m_];
            double temp = sum / jacobian_[in + in * m_];
            for (dsize i = in; i < m_; i++)
              jacobian_[i + k * m_] -= temp * jacobian_[i + in * m_];
            if (Rdiag[k] != 0.0) {
              temp = jacobian_[in + k * m_] / Rdiag[k];
              Rdiag[k] *= sqrt( std::max( 0.0, 1.0 - temp * temp ) );
              temp = Rdiag[k] / work1[k];
              DBGPRINT("\t\tQRFAC TEST: 0.5 * %g^2 <= %g\n", temp, machine_epsilon);
              if (0.05 * (temp * temp) <= machine_epsilon) {
                DBGPRINT("\t\tTEST PASSED\n");
                Rdiag[k] = VecNorm( jacobian_.begin() + (in1 + k * m_), m_ - in - 1 );
                work1[k] = Rdiag[k];
              }
            }
            DBGPRINT("QRFAC Rdiag[%lu]= %g\n", k+1, Rdiag[k]);
          }
        }
      }
      Rdiag[in] = -ajnorm;
    }      
    // -------------------------------------------
    // | END QRFAC
    // -------------------------------------------

    PrintMatrix("QR", n_, m_, jacobian_);
    for (dsize in = 0; in != n_; in++)
      DBGPRINT("\tRdiag[%lu]= %12.6g    acnorm[%lu]= %12.6g    ipvt[%lu]= %i\n",
             in+1, Rdiag[in], in+1, JcolNorm_[in], in+1, Jpvt[in]+1);

    double xnorm = 0.0;
    if ( currentIt == 0 ) {
      // First iteration. Scale according to the norms of the columns of
      // the initial Jacobian.
      for (dsize in = 0; in != n_; in++) {
        if ( JcolNorm_[in] == 0.0 )
          diag[in] = 1.0;
        else
          diag[in] = JcolNorm_[in];
      }
      PrintVector("diag", diag);
      // Calculate norm of scaled params and init step bound delta
      for (dsize in = 0; in != n_; in++)
        work1[in] = diag[in] * Params_[in];
      xnorm = VecNorm( work1 );
      delta = factor * xnorm;
      if (delta == 0.0)
        delta = factor;
      DBGPRINT("Delta= %g\n", delta);
    }

    // Form Qt * residual and store in Qt_r. Only first n components of
    // Qt*r are needed.
    Darray Qt_r( n_, 0.0 );
    newResidual = residual_;
    for (dsize in = 0; in != n_; in++) {
      double matElt = jacobian_[in + in * m_];
      DBGPRINT("DEBUG: Element %lu = %g\n", in+1, matElt);
      if (matElt != 0.0) {
        double sum = 0.0;
        for (dsize i = in; i < m_; i++)
          sum += jacobian_[i + in *  m_] * newResidual[i];
        double temp = -sum / matElt;
        for (dsize i = in; i < m_; i++)
          newResidual[i] += jacobian_[i + in * m_] * temp;
      }
      jacobian_[in + in * m_] = Rdiag[in]; 
      DBGPRINT("\tRdiag[%lu]= %g\n", in+1, jacobian_[in + in * m_]);
      Qt_r[in] = newResidual[in];
      DBGPRINT("DEBUG: qtf[%lu]= %g\n", in+1, Qt_r[in]);
    }

    // Compute the norm of the scaled gradient.
    double gnorm = 0.0;
    if ( rnorm > 0.0 ) {
      for (dsize in = 0; in != n_; in++) {
        int l = Jpvt[ in ];
        if (JcolNorm_[in] != 0.0) {
          double sum = 0.0;
          for (dsize i2 = 0; i2 <= in; i2++) {
            sum += jacobian_[i2 + in * m_] * (Qt_r[i2] / rnorm);
            DBGPRINT("DEBUG: jacobian[%lu, %lu]= %g\n", i2+1, in+1, jacobian_[i2 + in * m_]);
          }
          // Determine max
          double d1 = fabs( sum / JcolNorm_[l] );
          gnorm = std::max( gnorm, d1 );
        }
      }
    }
    DBGPRINT("gnorm= %g\n", gnorm);
    // Test for convergence of gradient norm.
    if (gnorm <= gtol) {
      DBGPRINT("Gradient norm %g is less than gtol %g, iteration %i\n",
             gnorm, gtol, currentIt+1);
      info = 4;
      break;
    }

    // Rescale if necessary
    for (dsize in = 0; in != n_; in++)
      diag[in] = std::max( diag[in], JcolNorm_[in] );
    PrintVector("RescaleDiag", diag);
    double Ratio = 0.0;
    while (Ratio < 0.0001) {
      // -----------------------------------------
      // | LMPAR BEGIN
      // -----------------------------------------
      // NOTE: Adapted from lmpar_ in lmdif.c from Grace 5.1.22
      DBGPRINT("\nLMPAR ITERATION %i\n", currentIt+1);
      // Determine the Levenberg-Marquardt parameter.
      Darray Xvec(n_, 0.0);  // Will contain solution to A*x=b, sqrt(par)*D*x=0
      // Compute and store in Xvec the Gauss-Newton direction.
      // If the Jacobian is rank-deficient, obtain a least-squares solution.
      dsize rank = n_;
      for (dsize in = 0; in != n_; in++) {
        work1[in] = Qt_r[in];
        DBGPRINT("jac[%lu,%lu]= %12.6g  rank= %4lu", in+1, in+1, jacobian_[in + in * m_], rank);
        if (jacobian_[in + in * m_] == 0.0 && rank == n_)
          rank = in;
        if (rank < n_)
          work1[in] = 0.0;
        DBGPRINT("   wa1= %12.6g\n", work1[in]);
      }
      DBGPRINT("Final rank= %lu\n", rank);
      // Subtract 1 from rank to use as an index.
      dsize rm1 = rank - 1;
      if (rm1 >= 0) {
        for (dsize k = 0; k <= rm1; k++) {
          dsize in = rm1 - k;
          work1[in] /= jacobian_[in + in * m_];
          DBGPRINT("wa1[%lu] /= %g\n", in+1, jacobian_[in + in * m_]);
          long int in1 = (long int)in - 1;
          if (in1 > -1) {
            double temp = work1[in];
            for (long int i = 0; i <= in1; i++) {
              work1[i] -= jacobian_[i + in * m_] * temp;
              DBGPRINT("  wa1[%li] -= %g\n", i+1, jacobian_[i + in * m_]);
            }
          }
        }
      }
      PrintVector("work1", work1);
      for (dsize in = 0; in != n_; in++)
        Xvec[ Jpvt[in] ] = work1[in];
  
      // Evaluate the function at the origin, and test for acceptance of
      // the Gauss-Newton direction.
      for (dsize in = 0; in != n_; in++) {
        DBGPRINT("work2[%lu] = %g * %g\n", in+1, diag[in], Xvec[in]);
        work2[in] = diag[in] * Xvec[in];
      }
      double dxnorm = VecNorm( work2 );
      DBGPRINT("dxnorm= %g\n", dxnorm);
      double fp = dxnorm - delta;
      // Initialize counter for searching for LM parameter
      int lmIterations = 0;
      if (fp > 0.1 * delta) {
        // If the Jacobian is not rank deficient, the Newton step provdes a lower
        // bound, parl, for the zero of the function. Otherwise set this bound to
        // zero
        double parl = 0.0;
        if ( rank >= n_ ) {
          for (dsize in = 0; in != n_; in++) {
            int idx = Jpvt[in];
            work1[in] = diag[idx] * (work2[idx] / dxnorm);
          }
          for (dsize in = 0; in != n_; in++) {
            double sum = 0.0;
            long int in1 = (long int)in - 1;
            if (in1 > -1) {
              for (long int i = 0; i <= in1; i++)
                sum += jacobian_[i + in * m_] * work1[i];
            }
            work1[in] = (work1[in] - sum) / jacobian_[in + in * m_];
          }
          double temp = VecNorm(work1);
          parl = fp / delta / temp / temp;
        }
        DBGPRINT("parl= %g\n",parl);
  
        // Calculate an upper bound, paru, for the zero of the function
        for (dsize in = 0; in != n_; in++) {
          double sum = 0.0;
          for (dsize i = 0; i <= in; i++)
            sum += jacobian_[i + in * m_] * Qt_r[i];
          work1[in] = sum / diag[ Jpvt[in] ];
          DBGPRINT("paru work1[%lu]= %g\n", in+1, work1[in]);
        }
        double w1norm = VecNorm( work1 );
        double paru = w1norm / delta;
        DBGPRINT("paru = %g = %g / %g\n", paru, w1norm, delta);
        if (paru == 0.0)
          paru = dwarf / std::min( delta, 0.1 );
        DBGPRINT("paru= %g\n", paru);
        
        // If the current L-M parameter lies outside of the interval (parl,paru),
        // set par to the closer endpoint.
        LM_par = std::max( LM_par, parl );
        LM_par = std::min( LM_par, paru );
        if (LM_par == 0.0)
          LM_par = w1norm / dxnorm;
  
        // Iteration start.
        bool lmLoop = true;
        while (lmLoop) {
          ++lmIterations;
          DBGPRINT("\t[ lmLoop %i  parl= %g  paru= %g   LM_par= %g ]\n",
                 lmIterations, parl, paru, LM_par);
          // Evaluate the function at the current value of LM_par
          if (LM_par == 0.0)
            LM_par = std::max( dwarf, 0.001 * paru );
          double temp = sqrt( LM_par );
          for (dsize in = 0; in != n_; in++) {
            work1[in] = temp * diag[in];
            DBGPRINT("\tLMPAR work1[%lu]= %g = %g * %g\n", in+1, work1[in], temp, diag[in]);
          }
  
          // -----------------------------------------------
          // | BEGIN QRSOLV
          // -----------------------------------------------
          // NOTE: Adapted from qrsolv_ in lmdif.c from Grace 5.1.22 
          DBGPRINT("\nQRSOLV ITERATION %i\n", lmIterations);
          // This array of length n will hold the diagonal elements of the
          // upper triangular matrix s.
          Darray Sdiag(n_, 0.0);

          // Copy R and Qt_r to preserve input and initialize s.
          for (dsize in = 0; in != n_; in++) {
            for (dsize i = in; i != n_; i++)
              jacobian_[i + in * m_] = jacobian_[in + i * m_];
            Xvec[in] = jacobian_[in + in * m_];
            work2[in] = Qt_r[in];
          }

          // Eliminate the diagonal matrix D using a Givens rotation
          for (dsize in = 0; in != n_; in++) {
            // Prepare the row of D to be eliminated, locating the diagonal 
            // element using p from the QR factorization.
            double diagL = work1[ Jpvt[in] ];
            DBGPRINT("diag[%i] = %g\n", Jpvt[in]+1, diagL);
            if (diagL != 0.0) {
              for (dsize k = in; k < n_; k++)
                Sdiag[k] = 0.0;
              Sdiag[in] = diagL;
              // The transformations to eliminate the row of D modify only
              // a single element of Qt_r beyond the first n, which is 
              // initially zero.
              double qtbpj = 0.0;
              for (dsize k = in; k < n_; k++) {
                // Determine a Givens rotation which eliminates the appropriate
                // element in the current row of D.
                if (Sdiag[k] != 0.0) {
                  double d1 = jacobian_[k + k * m_];
                  double cos, sin;
                  if (fabs(d1) < fabs(Sdiag[k])) {
                    double cotan = jacobian_[k + k * m_] / Sdiag[k];
                    sin = 0.5 / sqrt(0.25 + 0.25 * (cotan * cotan));
                    cos = sin * cotan;
                  } else {
                    double tan = Sdiag[k] / jacobian_[k + k *m_];
                    cos = 0.5 / sqrt(0.25 + 0.25 * (tan * tan));
                    sin = cos * tan;
                  }
                  DBGPRINT("DBG QRsolv: %lu, %lu: cos= %12.6g  sin= %12.6g\n",
                         in+1, k+1, cos, sin);
                  // Compute the modified diagonal element of R and the 
                  // modified element of (Qt_r, 0)
                  jacobian_[k + k * m_] = cos * jacobian_[k + k * m_] + sin * Sdiag[k];
                  double temp = cos * work2[k] + sin * qtbpj; 
                  qtbpj = -sin * work2[k] + cos * qtbpj;
                  work2[k] = temp;
                  DBGPRINT("\twork2[%lu]= %g\n", k+1, work2[k]);
                  // Accumulate the transformation in the row of S
                  dsize kp1 = k + 1;
                  if (kp1 <= n_) {
                    for (dsize i = kp1; i < n_; i++) {
                      temp = cos * jacobian_[i + k * m_] + sin * Sdiag[i];
                      Sdiag[i] = -sin * jacobian_[i + k * m_] + cos * Sdiag[i];
                      jacobian_[i + k * m_] = temp;
                    }
                  }
                }
              }
            }
            // Store the diagonal element of s and restore the corresponding
            // diagonal element of r
            Sdiag[in] = jacobian_[in + in * m_];
            DBGPRINT("QRsolv sdiag[%lu]= %g\n", in+1, Sdiag[in]);
            jacobian_[in + in * m_] = Xvec[in];
          }
  
          // Solve the triangular system for z. if the system is singular,
          // then obtain a least-squares solution.
          dsize nsing = n_;
          for (dsize in = 0; in != n_; in++) {
            if (Sdiag[in] == 0.0 && nsing == n_)
              nsing = in;
            if (nsing < n_)
              work2[in] = 0.0;
          }
          DBGPRINT("nsing= %lu\n", nsing);
          // Subtract 1 from nsing to use as an index.
          --nsing; 
          if (nsing >= 0) {
            for (dsize k = 0; k <= nsing; k++) {
              dsize in = nsing - k;
              double sum = 0.0;
              dsize in1 = in + 1;
              if (in1 <= nsing) {
                for (dsize i = in1; i <= nsing; i++)
                  sum += jacobian_[i + in * m_] * work2[i];
              }
              work2[in] = (work2[in] - sum) / Sdiag[in];
            }
          }

          // Permute the components of z back to components of x.
          for (dsize in = 0; in != n_; in++)
            Xvec[ Jpvt[in] ] = work2[in];

          PrintVector("QRsolv Xvec", Xvec);
          PrintVector("sdiag", Sdiag);
          // -----------------------------------------------
          // | END QRSOLV
          // -----------------------------------------------

          for (dsize in = 0; in != n_; in++) {
            work2[in] = diag[in] * Xvec[in];
            DBGPRINT("\twork2[%lu]= %g = %g * %g\n", in+1, work2[in], diag[in], Xvec[in]);
          }
  
          dxnorm = VecNorm( work2 );
          temp = fp;
          fp = dxnorm - delta;
          DBGPRINT("DBG LMPAR: fp = %g = %g - %g\n", fp, dxnorm, delta); 
          // If the function is small enough, accept the current value of LM_par.
          // Also test for the exceptional cases where parl is zero of the number
          // of iterations has reached 10.
          if (fabs(fp) <= 0.1 * delta ||
              (parl == 0.0 && fp <= temp && temp < 0.0) ||
              lmIterations == 10)
          {
            lmLoop = false;
            break;
          }
  
          // Compute the Newton correction.
          for (dsize in = 0; in != n_; in++) {
            int idx = Jpvt[ in ];
            work1[in] = diag[idx] * (work2[idx] / dxnorm);
          }
          for (dsize in = 0; in != n_; in++) {
            work1[in] /= Sdiag[in];
            temp = work1[in];
            dsize in1 = in + 1;
            if (in1 <= n_) {
              for (dsize i = in1; i < n_; i++)
                work1[i] -= jacobian_[i + in * m_] * temp;
            }
          }
          temp = VecNorm( work1 );
          double parc = fp / delta / temp / temp;
          DBGPRINT("parc= %g\n", parc);
          
          // Depending on the sign of the function, update parl or paru
          if (fp > 0.0)
            parl = std::max( parl, LM_par );
          if (fp < 0.0)
            paru = std::min( paru, LM_par );

          // Compute an improved estimate for the parameter
          LM_par = std::max( parl, LM_par + parc );
        }
      }
      DBGPRINT("END LMPAR iter= %i  LM_par= %g\n", lmIterations, LM_par); 
      if (lmIterations == 0)
        LM_par = 0.0;
      // -------------------------------------------
      // | LMPAR END
      // -------------------------------------------
      DBGPRINT("DEBUG: LM_par is %g\n", LM_par);
      // Store the direction Xvec and Param + Xvec. Calculate norm of Param.
      PrintVector("DEBUG: hvec", Xvec);
      for (dsize in = 0; in != n_; in++) {
        Xvec[in] = -Xvec[in];
        work1[in] = Params_[in] + Xvec[in];
        work2[in] = diag[in] * Xvec[in];
      }
      double pnorm = VecNorm( work2 );
      DBGPRINT("pnorm= %g\n", pnorm);

      // On first iteration, adjust initial step bound
      if (currentIt == 0)
        delta = std::min( delta, pnorm );
      DBGPRINT("Delta is now %g\n", delta);

      // Evaluate function at Param + Xvec and calculate its norm
      EvaluateFxn( Xvals_, Yvals_, work1, newResidual );
      ++nfev;
      double rnorm1 = VecNorm( newResidual );
      DBGPRINT("rnorm1= %g\n", rnorm1);

      // Compute the scaled actual reduction
      double actual_reduction = -1.0;
      if ( 0.1 * rnorm1 < rnorm) {
        double d1 = rnorm1 / rnorm;
        actual_reduction = 1.0 - d1 * d1;
      }
      DBGPRINT("actualReduction= %g\n", actual_reduction);

      // Compute the scaled predicted reduction and the scaled directional
      // derivative.
      for (dsize in = 0; in != n_; in++)
      {
        work2[in] = 0.0;
        double temp = Xvec[ Jpvt[in] ];
        for (dsize i = 0; i <= in; i++)
          work2[i] += jacobian_[i + in * m_] * temp;
      }
      double temp1 = VecNorm( work2 ) / rnorm;
      double temp2 = sqrt(LM_par) * pnorm / rnorm;
      DBGPRINT("temp1= %g    temp2= %g\n", temp1, temp2);

      double predicted_reduction = temp1 * temp1 + temp2 * temp2 / 0.5;
      double dirder = -(temp1 * temp1 + temp2 * temp2);
      DBGPRINT("predictedReduction = %12.6g    dirder= %12.6g\n", predicted_reduction, dirder);

      // Compute the ratio of the actial to the predicted reduction.
      if (predicted_reduction != 0.0)
        Ratio = actual_reduction / predicted_reduction;
      DBGPRINT("ratio= %g\n", Ratio);

      // Update the step bound
      if (Ratio <= 0.25) {
        DBGPRINT("Ratio <= 0.25\n");
        double temp;
        if (actual_reduction < 0.0)
          temp = 0.5 * dirder / (dirder + 0.5 * actual_reduction);
        else
          temp = 0.5;
        if ( 0.1 * rnorm1 >= rnorm || temp < 0.1 )
          temp = 0.1;
        DBGPRINT("delta = %g * min( %g, %g )\n", temp, delta, pnorm / 0.1); 
        delta = temp * std::min( delta, pnorm / 0.1 );
        LM_par /= temp;
      } else {
        if (LM_par == 0.0 || Ratio >= 0.75) {
          DBGPRINT("Ratio > 0.25 and (LMpar is zero or Ratio >= 0.75\n");
          delta = pnorm / 0.5;
          LM_par *= 0.5;
        }
      }
      DBGPRINT("LM_par= %g    Delta= %g\n", LM_par, delta);

      if (Ratio >= 0.0001) { 
        // Successful iteration. Update Param, residual, and their norms.
        for (dsize in = 0; in != n_; in++) {
          Params_[in] = work1[in];
          work1[in] = diag[in] * Params_[in];
        }
        for (dsize im = 0; im < m_; im++)
          residual_[im] = newResidual[im];
        xnorm = VecNorm( work1 );
        rnorm = rnorm1;
      }
      
      // Tests for convergence
      if (fabs(actual_reduction) <= tolerance &&
          predicted_reduction <= tolerance &&
          0.5 * Ratio <= 1.0)
        info = 1;
      if (delta <= xtol * xnorm)
        info = 2;  
      if (fabs(actual_reduction) <= tolerance &&
          predicted_reduction <= tolerance &&
          0.5 * Ratio <= 1.0 &&
          info == 2)
        info = 3;
      if (info != 0) break;
                 
      // Tests for stringent tolerance
      if (nfev >= maxfev)
        info = 5;
      if (fabs(actual_reduction) <= machine_epsilon &&
          predicted_reduction <= machine_epsilon &&
          0.5 * Ratio <= 1.0)
        info = 6;
      if (delta <= machine_epsilon * xnorm)
        info = 7;
      if (gnorm <= machine_epsilon)
        info = 8;
      if (info != 0) break;

    } // END inner loop
    if (info != 0) break;
    currentIt++;
  }
  // Final parameters and Y at final parameters
  Params_to_Pvec(ParamVec, Params_);
  fxn_(Xvals_, ParamVec, finalY_);
# ifdef DBG_CURVEFIT
  DBGPRINT("%s\n", Message(info));
  DBGPRINT("Exiting with info value = %i\n", info);
  for (dsize in = 0; in != n_; in++)
    DBGPRINT("\tParams[%lu]= %g\n", in, ParamVec[in]);
# endif
  return info;
}

// CurveFit::CalcMeanStdev()
void CurveFit::CalcMeanStdev(Darray const& Vals, double& mean, double& stdev) const
{
  mean = 0.0;
  for (Darray::const_iterator val = Vals.begin(); val != Vals.end(); ++val)
    mean += *val;
  mean /= (double)Vals.size();
  stdev = 0.0;
  if (Vals.size() > 1) {
    for (Darray::const_iterator val = Vals.begin(); val != Vals.end(); ++val)
    {
      double diff = mean - *val;
      stdev += (diff * diff);
    }
    stdev = sqrt( stdev / (double)(Vals.size() - 1) );
  }
}

// CurveFit::Statistics()
/** For final Y values vs given (presumably initial) Y values, calculate
  * correlation coefficient, chi-squared, Theil's U, and RMS percent 
  * error.
  */
int CurveFit::Statistics(Darray const& Yvals,
                         double& corr, double& ChiSq, 
                         double& TheilU, double& rms_percent_error) const
{
  if (finalY_.empty() || Yvals.size() != finalY_.size())
    return 9;
  int err_val = 0;
  unsigned int Nvals = Yvals.size();
  // Correlation coefficient
  corr = 0.0;
  if (Nvals > 1) {
    double mean0, mean1, stdev0, stdev1;
    CalcMeanStdev(finalY_, mean0, stdev0);
    CalcMeanStdev(Yvals,   mean1, stdev1);
    if (stdev0 > 0.0 && stdev1 > 0.0) {
      for (unsigned int n = 0; n != Nvals; n++)
        corr += (finalY_[n] - mean0) * (Yvals[n] - mean1);
      corr /= ((double)(Nvals - 1) * stdev0 * stdev1);
    }
  }
  ChiSq = 0.0;
  double Y2 = 0.0;
  bool setHasZero = false;
  for (unsigned int n = 0; n != Nvals; n++) {
    double diff = finalY_[n] - Yvals[n];
    ChiSq += (diff * diff);
    Y2 += (Yvals[n] * Yvals[n]);
    if (Yvals[n] == 0.0)
      setHasZero = true;
  }
  TheilU = sqrt(ChiSq / Y2);
  rms_percent_error = 0.0;
  if (!setHasZero) {
    for (unsigned int n = 0; n != Nvals; n++) {
      double diff = finalY_[n] - Yvals[n];
      rms_percent_error += (diff * diff) / (Yvals[n] * Yvals[n]);
    }
    rms_percent_error = sqrt( rms_percent_error / (double)Nvals );
  } else
    err_val = 10;
  return err_val;
}
