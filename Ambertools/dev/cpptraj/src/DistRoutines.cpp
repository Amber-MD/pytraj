#include <cmath> // sqrt
#include "DistRoutines.h"

/** \param a1 First set of XYZ coordinates.
  * \param a2 Second set of XYZ coordinates.
  * \param ucell Unit cell vectors.
  * \param recip Fractional cell vectors.
  * \return the shortest imaged distance^2 between the coordinates.
  */
double DIST2_ImageNonOrtho(Vec3 const& a1, Vec3 const& a2, 
                           Matrix_3x3 const& ucell, Matrix_3x3 const& recip) 
{ 
  int ixyz[3];
  return DIST2_ImageNonOrthoRecip(recip * a2, recip * a1, -1.0, ixyz, ucell);
}

/** \param a1 First set of XYZ coordinates.
  * \param a2 Second set of XYZ coordinates.
  * \param ucell Unit cell vectors.
  * \param recip Fractional cell vectors.
  * \return the shortest vector from a1 to a2. 
  */
Vec3 MinImagedVec(Vec3 const& a1, Vec3 const& a2,
                  Matrix_3x3 const& ucell, Matrix_3x3 const& recip)
{
  Vec3 f1 = recip * a1;
  Vec3 f2 = recip * a2;
//  mprintf("DEBUG: a1= %g %g %g, f1= %g %g %g\n", a1[0], a1[1], a1[2], f1[0], f1[1], f1[2]);
//  mprintf("DEBUG: a2= %g %g %g, f2= %g %g %g\n", a2[0], a2[1], a2[2], f2[0], f2[1], f2[2]);
  for (unsigned int i = 0; i < 3; i++) {
    f1[i] = f1[i] - floor(f1[i]);
    f2[i] = f2[i] - floor(f2[i]);
  }
  // Self
  Vec3 c1 = ucell.TransposeMult( f1 );
  Vec3 c2 = ucell.TransposeMult( f2 );
  Vec3 minVec = c2 - c1;
  double minDist2 = minVec.Magnitude2();
  // Images
  for (int ix = -1; ix < 2; ix++) {
    for (int iy = -1; iy < 2; iy++) {
      for (int iz = -1; iz < 2; iz++) {
        if (ix != 0 || iy != 0 || iz != 0) { // Ignore a2 self
          Vec3 ixyz(ix, iy, iz);
          c2 = ucell.TransposeMult(f2 + ixyz); // a2 image back in Cartesian space
          Vec3 dxyz = c2 - c1;
          double dist2 = dxyz.Magnitude2();
          if (dist2 < minDist2) {
            minDist2 = dist2;
            minVec = dxyz;
          }
        }
      }
    }
  }
  return minVec;
}

/** NON-ORTHORHOMBIC CASE: find shortest distance in periodic reference
  * This is a brute force check requiring up to 26 distance evaluations.
  * It has been adapted to be smarter by returning the first distance that
  * is shorter than the minimum possible distance between images.
  * \param f1 First set of fractional XYZ coordinates.
  * \param f2 Second set of fractional XYZ coordinates.
  * \param ixyz Will be set with integer coefficients describing closest reflection.
  * \param ucell Unit cell vectors.
  * \return the shortest imaged distance^2 between the coordinates.
  */
double DIST2_ImageNonOrthoRecip(Vec3 const& f1, Vec3 const& f2, double minIn, 
                                int* ixyz, Matrix_3x3 const& ucell) 
{ 
  //double closest2
  // The floor() calls serve to bring each point back in the main unit cell.
  double fx = f1[0] - floor(f1[0]);
  double fy = f1[1] - floor(f1[1]);
  double fz = f1[2] - floor(f1[2]); 
  double f2x = f2[0] - floor(f2[0]);
  double f2y = f2[1] - floor(f2[1]);
  double f2z = f2[2] - floor(f2[2]);
  // f2 back in Cartesian space
  double X_factor = (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
  double Y_factor = (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
  double Z_factor = (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);
  // Precompute some factors
  double fxm1 = fx - 1.0;
  double fxp1 = fx + 1.0;
  double fym1 = fy - 1.0;
  double fyp1 = fy + 1.0;
  double fzm1 = fz - 1.0;
  double fzp1 = fz + 1.0;

  double fxm1u0 = fxm1 * ucell[0];
  double fxu0   = fx   * ucell[0];
  double fxp1u0 = fxp1 * ucell[0];
  double fxm1u1 = fxm1 * ucell[1];
  double fxu1   = fx   * ucell[1];
  double fxp1u1 = fxp1 * ucell[1];
  double fxm1u2 = fxm1 * ucell[2];
  double fxu2   = fx   * ucell[2];
  double fxp1u2 = fxp1 * ucell[2];

  double fym1u3 = fym1 * ucell[3];
  double fyu3   = fy   * ucell[3];
  double fyp1u3 = fyp1 * ucell[3];
  double fym1u4 = fym1 * ucell[4];
  double fyu4   = fy   * ucell[4];
  double fyp1u4 = fyp1 * ucell[4];
  double fym1u5 = fym1 * ucell[5];
  double fyu5   = fy   * ucell[5];
  double fyp1u5 = fyp1 * ucell[5];

  double fzm1u6 = fzm1 * ucell[6];
  double fzu6   = fz   * ucell[6];
  double fzp1u6 = fzp1 * ucell[6];
  double fzm1u7 = fzm1 * ucell[7];
  double fzu7   = fz   * ucell[7];
  double fzp1u7 = fzp1 * ucell[7];
  double fzm1u8 = fzm1 * ucell[8];
  double fzu8   = fz   * ucell[8];
  double fzp1u8 = fzp1 * ucell[8];

  // Calc ix iy iz = 0 case
  double x = (fxu0 + fyu3 + fzu6) - X_factor;
  double y = (fxu1 + fyu4 + fzu7) - Y_factor;
  double z = (fxu2 + fyu5 + fzu8) - Z_factor;
  // DEBUG
  //mprintf("DEBUG: a2: %g %g %g\n",(fxu0 + fyu3 + fzu6), (fxu1 + fyu4 + fzu7), (fxu2 + fyu5 + fzu8));
  //mprintf("DEBUG: a1: %g %g %g\n", X_factor, Y_factor, Z_factor);
  double min = (x*x) + (y*y) + (z*z);

  if (minIn > 0.0 && minIn < min) min = minIn;

  ixyz[0] = 0;
  ixyz[1] = 0;
  ixyz[2] = 0;

  // -1 -1 -1
  x = (fxm1u0 + fym1u3 + fzm1u6) - X_factor;
  y = (fxm1u1 + fym1u4 + fzm1u7) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzm1u8) - Z_factor;
  double D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] = -1; }
  // -1 -1  0
  x = (fxm1u0 + fym1u3 + fzu6  ) - X_factor;
  y = (fxm1u1 + fym1u4 + fzu7  ) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] =  0; }
  // -1 -1 +1
  x = (fxm1u0 + fym1u3 + fzp1u6) - X_factor;
  y = (fxm1u1 + fym1u4 + fzp1u7) - Y_factor;
  z = (fxm1u2 + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] = -1; ixyz[2] =  1; }
  // -1  0 -1
  x = (fxm1u0 + fyu3   + fzm1u6) - X_factor;
  y = (fxm1u1 + fyu4   + fzm1u7) - Y_factor;
  z = (fxm1u2 + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] = -1; }
  // -1  0  0
  x = (fxm1u0 + fyu3   + fzu6  ) - X_factor;
  y = (fxm1u1 + fyu4   + fzu7  ) - Y_factor;
  z = (fxm1u2 + fyu5   + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] =  0; }
  // -1  0 +1
  x = (fxm1u0 + fyu3   + fzp1u6) - X_factor;
  y = (fxm1u1 + fyu4   + fzp1u7) - Y_factor;
  z = (fxm1u2 + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  0; ixyz[2] =  1; }
  // -1 +1 -1
  x = (fxm1u0 + fyp1u3 + fzm1u6) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] = -1; }
  // -1 +1  0
  x = (fxm1u0 + fyp1u3 + fzu6  ) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] =  0; }
  // -1 +1 +1
  x = (fxm1u0 + fyp1u3 + fzp1u6) - X_factor;
  y = (fxm1u1 + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxm1u2 + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] = -1; ixyz[1] =  1; ixyz[2] =  1; }

  //  0 -1 -1
  x = (fxu0   + fym1u3 + fzm1u6) - X_factor;
  y = (fxu1   + fym1u4 + fzm1u7) - Y_factor;
  z = (fxu2   + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] = -1; }
  //  0 -1  0
  x = (fxu0   + fym1u3 + fzu6  ) - X_factor;
  y = (fxu1   + fym1u4 + fzu7  ) - Y_factor;
  z = (fxu2   + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] =  0; }
  //  0 -1 +1
  x = (fxu0   + fym1u3 + fzp1u6) - X_factor;
  y = (fxu1   + fym1u4 + fzp1u7) - Y_factor;
  z = (fxu2   + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] = -1; ixyz[2] =  1; }
  //  0  0 -1
  x = (fxu0   + fyu3   + fzm1u6) - X_factor;
  y = (fxu1   + fyu4   + fzm1u7) - Y_factor;
  z = (fxu2   + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  0; ixyz[2] = -1; }
  //  0  0  0
  //  0  0 +1
  x = (fxu0   + fyu3   + fzp1u6) - X_factor;
  y = (fxu1   + fyu4   + fzp1u7) - Y_factor;
  z = (fxu2   + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  0; ixyz[2] =  1; }
  //  0 +1 -1
  x = (fxu0   + fyp1u3 + fzm1u6) - X_factor;
  y = (fxu1   + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxu2   + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] = -1; }
  //  0 +1  0
  x = (fxu0   + fyp1u3 + fzu6  ) - X_factor;
  y = (fxu1   + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxu2   + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] =  0; }
  //  0 +1 +1
  x = (fxu0   + fyp1u3 + fzp1u6) - X_factor;
  y = (fxu1   + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxu2   + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  0; ixyz[1] =  1; ixyz[2] =  1; }

  // +1 -1 -1
  x = (fxp1u0 + fym1u3 + fzm1u6) - X_factor;
  y = (fxp1u1 + fym1u4 + fzm1u7) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] = -1; }
  // +1 -1  0
  x = (fxp1u0 + fym1u3 + fzu6  ) - X_factor;
  y = (fxp1u1 + fym1u4 + fzu7  ) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] =  0; }
  // +1 -1 +1
  x = (fxp1u0 + fym1u3 + fzp1u6) - X_factor;
  y = (fxp1u1 + fym1u4 + fzp1u7) - Y_factor;
  z = (fxp1u2 + fym1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] = -1; ixyz[2] =  1; }
  // +1  0 -1
  x = (fxp1u0 + fyu3   + fzm1u6) - X_factor;
  y = (fxp1u1 + fyu4   + fzm1u7) - Y_factor;
  z = (fxp1u2 + fyu5   + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] = -1; }
  // +1  0  0
  x = (fxp1u0 + fyu3   + fzu6  ) - X_factor;
  y = (fxp1u1 + fyu4   + fzu7  ) - Y_factor;
  z = (fxp1u2 + fyu5   + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] =  0; }
  // +1  0 +1
  x = (fxp1u0 + fyu3   + fzp1u6) - X_factor;
  y = (fxp1u1 + fyu4   + fzp1u7) - Y_factor;
  z = (fxp1u2 + fyu5   + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  0; ixyz[2] =  1; }
  // +1 +1 -1
  x = (fxp1u0 + fyp1u3 + fzm1u6) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzm1u7) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzm1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] = -1; }
  // +1 +1  0
  x = (fxp1u0 + fyp1u3 + fzu6  ) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzu7  ) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzu8  ) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] =  0; }
  // +1 +1 +1
  x = (fxp1u0 + fyp1u3 + fzp1u6) - X_factor;
  y = (fxp1u1 + fyp1u4 + fzp1u7) - Y_factor;
  z = (fxp1u2 + fyp1u5 + fzp1u8) - Z_factor;
  D = (x*x) + (y*y) + (z*z);
  if (D < min) { min = D; ixyz[0] =  1; ixyz[1] =  1; ixyz[2] =  1; }

  //if (closest2 != 0.0 && min < closest2) return (min);
//  this->ClosestImage(a1, a2, ixyz);
//  fprintf(stdout,"DEBUG: Predict  = %2i %2i %2i\n",ixyz[0],ixyz[1],ixyz[2]);

//  ix = ixyz[0];
//  iy = ixyz[1];
//  iz = ixyz[2];

//D = sqrt(min);
//  fprintf(stdout,"DEBUG: MinDist  = %2i %2i %2i = %8.3f\n", ixmin, iymin, izmin, D);
//  printf("---------------------------------------------------------------\n");
  return(min);
}

// Frame::DIST2_ImageOrtho()
/** Return the minimum orthorhombic imaged distance^2 between coordinates a1 
  * and a2.
  */
double DIST2_ImageOrtho(Vec3 const& a1, Vec3 const& a2, Box const& box) {
  // If box lengths are zero no imaging possible
  if (box[0]==0.0 || box[1]==0.0 || box[2]==0.0) return -1.0;
  double x = a1[0] - a2[0];
  double y = a1[1] - a2[1];
  double z = a1[2] - a2[2];
  // Get rid of sign info
  if (x<0) x=-x;
  if (y<0) y=-y;
  if (z<0) z=-z;
  // Get rid of multiples of box lengths 
  while (x > box[0]) x = x - box[0];
  while (y > box[1]) y = y - box[1];
  while (z > box[2]) z = z - box[2];
  // Find shortest distance in periodic reference
  double D = box[0] - x;
  if (D < x) x = D;
  D = box[1] - y;
  if (D < y) y = D;  
  D = box[2] - z;
  if (D < z) z = D;

  return (x*x + y*y + z*z);
}

// Frame::DIST2_NoImage()
/** Return distance^2 between coordinates in a1 and a2.
  */
double DIST2_NoImage(const double* a1, const double* a2) {
  double x = a1[0] - a2[0];
  double y = a1[1] - a2[1];
  double z = a1[2] - a2[2];
  //double D = x*x + y*y + z*z;
  //fprintf(stdout,"Mask1=%8.3f %8.3f %8.3f Mask2=%8.3f %8.3f %8.3f D=%8.3f\n",
  //        a1[0],a1[1],a1[2],a2[0],a2[1],a2[2],D);
  return (x*x + y*y + z*z);
}

double DIST2_NoImage( Vec3 const& a1, Vec3 const& a2 ) {
  Vec3 vec = a1 - a2;
  return vec.Magnitude2();
}

double DIST_NoImage( Vec3 const& a1, Vec3 const& a2 ) {
  Vec3 vec = a1 - a2;
  return sqrt( vec.Magnitude2() );
}

double DIST2(const double* a1, const double* a2, ImagingType itype,
             Box const& box, Matrix_3x3 const& ucell, Matrix_3x3 const& recip)
{
  if (itype==NOIMAGE) 
    return DIST2_NoImage( a1, a2 );
  else if (itype==ORTHO) 
    return DIST2_ImageOrtho( a1, a2, box );
  else // NONORTHO
    return DIST2_ImageNonOrtho( a1, a2, ucell, recip );
}
