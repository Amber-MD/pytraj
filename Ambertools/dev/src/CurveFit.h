#ifndef INC_CURVEFIT_H
#define INC_CURVEFIT_H
#ifdef DBG_CURVEFIT
# include <cstdio> // DEBUG
#endif
#include <vector>
/// Used for non-linear curve fitting.
class CurveFit {
  public:
    typedef std::vector<double> Darray;
    /** Function prototype for Y = f(X):
      * X value, P[] parameters, dYdP[] first derivatives of Y w.r.t. P[]
      * \return Y
      */
    //typedef double (*FitFunctionType)(double, Darray const&, Darray&);
    /** Function prototype for Y =f(x):
      * Xvalue, P[] parameters.
      */
    //typedef double (*FitFunctionType)(double, Darray const&);
    /** Function prototype for Y = f(x, p)
      * Xvalues, P[] parameters, Y values
      */
    typedef int (*FitFunctionType)(Darray const&, Darray const&, Darray&);
    
    CurveFit();
#   ifdef DBG_CURVEFIT
    ~CurveFit(); // DEBUG
#   endif
    /// Perform Levenberg-Marquardt curve fit: fxn, x, y, p, tol, iter
    int LevenbergMarquardt(FitFunctionType, Darray const&, Darray const&, Darray&,
                           double, int);
    /// Perform Legenberg-Marquardt curve fit: fxn, x, y, p, wgt, tol, iter
    int LevenbergMarquardt(FitFunctionType, Darray const&, Darray const&, Darray&,
                           Darray const&, double, int);
    /// Perform Levenberg-Marquardt curve fit: fxn, x, y, p, bnd, lbnd, ubnd, wts, tol, iter
    int LevenbergMarquardt(FitFunctionType, Darray const&, Darray const&, Darray&,
                           std::vector<bool> const&, Darray const&, Darray const&, 
                           Darray const&, double, int);
    /// Calculate various statistics using final Y values and input Y values.
    int Statistics(Darray const&, double&, double&, double&, double&) const;
    /// \return Status message.
    static const char* Message(int);
    /// \return Error message if status is zero.
    const char* ErrorMessage() const { return errorMessage_; } 
    Darray const& FinalY() const { return finalY_; }
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<double>::size_type dsize;
    static const double machine_epsilon;
    /// Function params to internal params.
    void Pvec_to_Params(Darray&);
    /// Internal params to params for function evaluation.
    void Params_to_Pvec(Darray&, Darray const&) const;
    /// Calculate residual using given X, Y, and parameters.
    void EvaluateFxn(Darray const&, Darray const&, Darray const&, Darray&);
    /// Calculate Jacobian using forward-difference approximation
    void CalcJacobian_ForwardDiff(Darray const&, Darray const&, Darray&, Darray const&, Darray&);
    /// Calculate || m(i,...) || for row vector in matrix
    static double VecNorm(Darray::const_iterator const&, dsize);
    /// Calculate || v || for vector
    static inline double VecNorm(Darray const& vec) {
      return VecNorm( vec.begin(), vec.size() );
    }
    /// Print final parameters to STDOUT
    void PrintFinalParams(Darray const&) const;
    /// \return true if input coords/parameters have problems.
    bool ParametersHaveProblems(Darray const&, Darray const&, Darray const&);
    /// Calculate mean and standard deviation
    void CalcMeanStdev(Darray const&, double&, double&) const;
    // DEBUG
    inline void PrintMatrix(const char*, int, int, Darray const&) const;
    inline void PrintVector(const char*, Darray const&) const;
    void DBGPRINT(const char*, ...) const;

    FitFunctionType fxn_; ///< Function to fit to.
    dsize m_;             ///< Number of values (rows)
    dsize n_;             ///< Number of parameters (cols)
    Darray jacobian_; ///< Jacobian/R, stored in transpose (row-major)
    Darray Params_;   ///< Working copy of parameter vector.
    Darray fParms_;   ///< Parameters for function evaluation.
    Darray finalY_;     ///< fxn(x) with final parameters
    Darray Weights_;  ///< Residual weights
    std::vector<bool> hasBounds_; ///< True if parameter has bounds.
    Darray Ubound_;               ///< Parameter upper bound.
    Darray Lbound_;               ///< Parameter lower bound.
    const char* errorMessage_; ///< Set to error message when return status is 0.
    // DEBUG
#   ifdef DBG_CURVEFIT
    FILE* dbgfile_;
#   endif
};
#endif
