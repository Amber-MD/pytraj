#include "SimplexMin.h"
#include "CpptrajStdio.h"

// ---------- SIMPLEX MINIMIZER ------------------------------------------------
double SimplexMin::chi_squared(Darray const& Ysearch) { 
  double chisq = 0.0;

  fxn_(Xvals_, Ysearch, Ynew_);

  for (dsize i = 0; i < Nvals_; i++) {
    double diff = Yvals_[i] - Ynew_[i];
    chisq += (diff * diff);
  }

  return chisq;
}

// Amotry()
double SimplexMin::Amotry(Darray& psum, int ihi, double fac)
{
  Darray ptry(NP_);

  double fac1 = (1.0 - fac) / (double)NP_;
  double fac2 = fac1 - fac;
  for (dsize j = 0; j < NP_; j++)
  {
    ptry[j] = (psum[j] * fac1) - (Xsimplex_[ihi * NP_ + j] * fac2);
    //mprintf("\t\tAmotry: %6lu%10.5g\n",j,ptry[j]);
  }
  double ytry = chi_squared( ptry );
  if (ytry < Ysearch_[ihi]) {
    Ysearch_[ihi] = ytry;
    for (dsize j = 0; j < NP_; j++) {
      psum[j] = psum[j] - Xsimplex_[ihi * NP_ + j] + ptry[j];
      Xsimplex_[ihi * NP_ + j] = ptry[j];
      //mprintf("\t\tAmotryX: %6lu%6lu%10.5g\n",ihi+1,j+1,xsmplx[ihi][j]);
    }
  }
  return ytry;
}

// Amoeba()
/** Main driver for the simplex method */
int SimplexMin::Amoeba(int amoeba_itmax, double amoeba_ftol) {
  Darray psum(NP_); // TODO: Make class var?

  int iter = 0;
  bool loop1 = true;
  while (loop1) {
//    mprintf("Hit loop one %6i\n",iter);
    for (dsize n = 0; n < NP_; n++) {
      psum[n] = 0;
      for (dsize m = 0; m < NP1_; m++) { 
//        mprintf("Xsmplx %6lu%6lu%10.5g\n",m,n,Xsimplex_[m * NP_ + n]);
//        mprintf("Ysearch %6lu%10.5g\n",m,Ysearch_[m]);
        psum[n] += Xsimplex_[m * NP_ + n];
      }
    }
    bool loop2 = true;
    while (loop2) {
//      mprintf("Hit loop two %6i\n",iter);
//      for (dsize n = 0; n < NP_; n++)
//        mprintf("\tPsum %6lu%10.5g\n",n,psum[n]);
      dsize ilo = 0, ihi, inhi;
      if (Ysearch_[0] > Ysearch_[1]) {
        ihi = 0;
        inhi = 1;
      } else {
        ihi = 1;
        inhi = 0;
      }
      for (dsize i = 0; i < NP1_; i++) {
        if (Ysearch_[i] <= Ysearch_[ilo]) ilo = i;
        if (Ysearch_[i] > Ysearch_[ihi]) {
          inhi = ihi;
          ihi = i;
        } else if ( Ysearch_[i] > Ysearch_[inhi]) {
          if (i != ihi) inhi = i;
        }
      }  
//      mprintf("Yihi Yilo = %10.5g%10.5g\n",Ysearch_[ihi],Ysearch_[ilo]);
//      mprintf("\tYihi Yilo = %6lu%6lu\n",ihi+1,ilo+1);
      double abs_yhi_ylo = Ysearch_[ihi] - Ysearch_[ilo];
      if (abs_yhi_ylo < 0) abs_yhi_ylo = -abs_yhi_ylo;
      double abs_yhi = Ysearch_[ihi];
      if (abs_yhi < 0) abs_yhi = -abs_yhi;
      double abs_ylo = Ysearch_[ilo];
      if (abs_ylo < 0) abs_ylo = -abs_ylo;
//      mprintf("Abs(yihi - yilo)=%g, Abs(yihi)=%g, Abs(yilo)=%g\n",
//              abs_yhi_ylo,abs_yhi,abs_ylo);
      double rtol = 2.0 * (abs_yhi_ylo / (abs_yhi + abs_ylo));
      if (rtol < amoeba_ftol) {
        double swap = Ysearch_[0];
        Ysearch_[0] = Ysearch_[ilo];
        Ysearch_[ilo] = swap;
        for (dsize n = 0; n < NP_; n++) {
          swap = Xsimplex_[n]; // 0 * NP_ + n
          Xsimplex_[n] = Xsimplex_[ilo * NP_ + n];
          Xsimplex_[ilo * NP_ + n] = swap;
        }
        mprintf("Rtol %g is less than specified tolerance %g\n", rtol, amoeba_ftol);
        return iter;
      }
//      mprintf("\tIn amoeba, iter=%i, rtol=%15.6g\n",iter,rtol); 

      if (iter >= amoeba_itmax) {
        mprintf("Max iterations (%i) exceeded in amoeba.\n",amoeba_itmax);
        return iter;
      }
      iter += 2;

      double ytry = Amotry(psum, ihi, -1.0);
//      mprintf("\tYtry %6i%10.5g\n",iter,ytry);
      if (ytry <= Ysearch_[ilo]) { 
        ytry = Amotry( psum, ihi, 2.0);
//        mprintf("\tCase 1 %10.5g\n",ytry);
      } else if (ytry >= Ysearch_[inhi]) {
        double ysave = Ysearch_[ihi];
        ytry = Amotry(psum, ihi, 0.5);
//        mprintf("\tCase 2 %10.5g\n",ytry);
        if (ytry >= ysave) {
          for (dsize i=0; i < NP1_; i++) {
            if (i != ilo) {
              for (dsize j = 0; j < NP_; j++) {
                psum[j] = 0.5 * (Xsimplex_[i * NP_ + j] + Xsimplex_[ilo * NP_ + j]);
                //psum[j] = 0.5 * (xsmplx[i][j] + xsmplx[ilo][j]);
                Xsimplex_[i * NP_ + j] = psum[j];
                //xsmplx[i][j] = psum[j];
              }
              Ysearch_[i] = chi_squared(psum);
            }
          }
          iter += (int)NP_;
          // GO TO 1
          loop2 = false;
        }
      } else {
//        mprintf("\tCase 3\n");
        iter--;
      }
      // GO TO 2
    } // END LOOP 2
  } // END LOOP 1
  return iter;
}

// Average_vertices()
void SimplexMin::Average_vertices(Darray& xsearch) const
{
  for (dsize j=0; j < NP_; j++) {
    xsearch[j] = 0.0;
    for (dsize k = 0; k < NP1_; k++) 
      xsearch[j] += Xsimplex_[k * NP_ + j];
    xsearch[j] /= (double)NP1_;
  }
}

// SimplexMin::Simplex_min()
/** Main driver routine for Amoeba (downhill simplex) minimizer. In the 
  * simplex method, N+1 initial points (where N is the dimension of the 
  * search space) must be chosen. Initial points are stored in rows of
  * matrix Xsimplex_; one of these is the input Q_vector. The other  
  * initial solutions should be of the order of the characteristic 
  * "lengthscales" over which Q_vector varies. delqfracIn determines the 
  * size of variation for each of the components of Q; the sign of the 
  * variation is randomly chosen.
  */
int SimplexMin::Minimize(SimplexFunctionType fxnIn, Darray& Q_vector, 
                         DataSet* Xin, Darray const& YvalsIn, double delqfracIn,
                         int maxItIn, double ftolIn, int nsearchIn,
                         Random_Number& RNgen) // TODO: use internal rand
{
  double delqfrac = delqfracIn;
  Xvals_ = Xin;
  Yvals_ = YvalsIn;
  Ynew_ = YvalsIn;
  Nvals_ = Yvals_.size();
  fxn_ = fxnIn;

  NP_ = Q_vector.size();
  NP1_ = NP_ + 1;

  Ysearch_.assign(NP1_, 0.0);

  Xsimplex_.assign(NP1_ * NP_, 0.0);

  //int test_seed = -3001796; // For tensorfit_ comparison

  // Initial chi squared
  double chisq = chi_squared( Q_vector );
  mprintf("\tInitial chi-squared is %g\n", chisq);

  // Now execute the simplex search method with initial vertices,
  // xsimplx, and chi-squared values, ysearch.
  // We restart the minimization a pre-set number of times to avoid
  // an anomalous result.
  Darray xsearch = Q_vector;

  // BEGIN SIMPLEX LOOP
  for (int i = 0; i < nsearchIn; i++) {
    for (dsize j = 0; j < NP_; j++) 
      Xsimplex_[j] = xsearch[j]; // 0 * NP_ + j

    for (dsize j = 0; j < NP_; j++) {
      for (dsize k = 0; k < NP_; k++) {
        if (j==k) {
          double sgn = RNgen.rn_gen() - 0.5; // -0.5 <= sgn <= 0.5
          //sgn = random_(test_seed) - 0.5; // For tensorfit_ comparison
          if (sgn < 0)
            sgn = -1.0;
          else
            sgn = 1.0;
          Xsimplex_[(j+1) * NP_ + k] = Xsimplex_[k] * (1+(sgn*delqfrac));
          //xsmplx[j+1][k] = xsmplx[0][k] * (1+(sgn*delqfrac_));
        } else {
          Xsimplex_[(j+1) * NP_ + k] = Xsimplex_[k];
          //xsmplx[j+1][k] = xsmplx[0][k];
        }
      }
    }
//    mprintf("--------------------\n");
//    for (dsize j = 0; j < NP1_; j++)
//      for (dsize k = 0; k < NP_; k++)
//        mprintf("xsimplex[%lu,%lu]= %g\n", j, k, Xsimplex_[j * NP_ + k]);
//    for (dsize k = 0; k < NP_; k++)
//      mprintf("xsearch[%lu]= %g\n", k, xsearch[k]);
//    mprintf("--------------------\n");

    // As to amoeba, chi-squared must be evaluated for all
    // vertices in the initial simplex.
    for (dsize j=0; j < NP1_; j++) {
      for (dsize k=0; k < NP_; k++) {
        xsearch[k] = Xsimplex_[j * NP_ + k];
        //xsearch[k] = xsmplx[j][k];
      }
      Ysearch_[j] = chi_squared(xsearch);
//      mprintf("ysearch[%lu]= %g\n", j, Ysearch_[j]);
    }

    // Average the vertices and compute details of the average.
    Average_vertices( xsearch );
    chisq  = chi_squared( xsearch );

    mprintf("Input to amoeba - average at cycle %i\n",i+1);
    mprintf("    Initial chisq = %15.5g\n",chisq);

    int am_iter = Amoeba( maxItIn, ftolIn );
    mprintf("amoeba ran for %i iterations.\n", am_iter);
    // Put amoeba results into xsearch
    for (dsize j=0; j < NP1_; j++) {
      for (dsize k=0; k < NP_; k++) {
        //xsearch[k] = xsmplx[j][k];
        xsearch[k] = Xsimplex_[j * NP_ + k];
      }
      Ysearch_[j] = chi_squared(xsearch);
    }

    // Average the vertices and compute details of the average.
    Average_vertices( xsearch );
    final_chisq_ = chi_squared( xsearch );
    mprintf("Output from amoeba - average at cycle %i\n",i+1);
    mprintf("    Final chisq = %15.5g\n",final_chisq_);
   
    // cycle over main loop, but first reduce the size of delqfrac:
    delqfrac *= 0.750;
    mprintf("\tAmoeba: Setting delqfrac to %15.7g\n",delqfrac);
  }

  // Set q vector to the final average result from simpmin
  Q_vector = xsearch;

  return 0;
}
