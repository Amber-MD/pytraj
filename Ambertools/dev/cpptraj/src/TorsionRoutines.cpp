/*! \file TorsionRoutines.cpp
    \brief Routines used to calculate torsions and angles.
 */
#include <cmath>
#include "TorsionRoutines.h"
#include "Constants.h" // PI, TWOPI
#include "Vec3.h"

// Torsion()
/** Given 4 sets of XYZ coords, calculate the torsion (in radians) between the 
  * planes formed by a1-a2-a3 and a2-a3-a4.
  */
double Torsion(const double *a1, const double *a2, const double *a3, const double *a4) 
{
  double Lx = ((a2[1]-a1[1])*(a3[2]-a2[2])) - ((a2[2]-a1[2])*(a3[1]-a2[1])); 
  double Ly = ((a2[2]-a1[2])*(a3[0]-a2[0])) - ((a2[0]-a1[0])*(a3[2]-a2[2])); 
  double Lz = ((a2[0]-a1[0])*(a3[1]-a2[1])) - ((a2[1]-a1[1])*(a3[0]-a2[0]));

  double Rx = ((a4[1]-a3[1])*(a2[2]-a3[2])) - ((a4[2]-a3[2])*(a2[1]-a3[1])); 
  double Ry = ((a4[2]-a3[2])*(a2[0]-a3[0])) - ((a4[0]-a3[0])*(a2[2]-a3[2])); 
  double Rz = ((a4[0]-a3[0])*(a2[1]-a3[1])) - ((a4[1]-a3[1])*(a2[0]-a3[0]));

  double Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
  double Rnorm = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);

  double Sx = (Ly*Rz) - (Lz*Ry); 
  double Sy = (Lz*Rx) - (Lx*Rz); 
  double Sz = (Lx*Ry) - (Ly*Rx);

  double angle = (Lx*Rx + Ly*Ry + Lz*Rz) / (Lnorm * Rnorm);

  if ( angle > 1.0 ) angle = 1.0;
  if ( angle < -1.0 ) angle = -1.0;

  angle = acos( angle );

  if ( (Sx * (a3[0]-a2[0]) + Sy * (a3[1]-a2[1]) + Sz * (a3[2]-a2[2])) < 0 )
    angle = -angle;

  return angle;
}

/// Constant used in AS pucker calc
static const double pi_over_5 = Constants::PI / 5.0;

// Pucker_AS()
/** Return the pucker (in radians) of coords stored in a1-a5 based on 
  * Altona & Sundaralingam method.
  */
double Pucker_AS(const double* a1, const double* a2, const double* a3, 
                 const double* a4, const double* a5, double& amp) 
{
  double pucker;
  double v1, v2, v3, v4, v5, a, b;

  pucker = 0.0;
  amp = 0.0;

  v4 = Torsion(a4, a5, a1, a2);
  v5 = Torsion(a5, a1, a2, a3);
  v1 = Torsion(a1, a2, a3, a4);
  v2 = Torsion(a2, a3, a4, a5);
  v3 = Torsion(a3, a4, a5, a1);

  a = (v1*cos(0.0) +
       v2*cos( 4.0*pi_over_5) +
       v3*cos( 8.0*pi_over_5) +
       v4*cos(12.0*pi_over_5) +
       v5*cos(16.0*pi_over_5))*0.4;

  b = (v1*sin(0.0) +
       v2*sin( 4.0*pi_over_5) +
       v3*sin( 8.0*pi_over_5) +
       v4*sin(12.0*pi_over_5) +
       v5*sin(16.0*pi_over_5))*-0.4;

  amp = sqrt(a*a + b*b);

  if (amp != 0.0)
    pucker = atan2(b,a);
  if (pucker < 0) pucker += Constants::TWOPI;

  return pucker;
}

// Pucker_CP()
/** Return the pucker (in radians) of coords in a1-a5 based on method of
  * Cremer & Pople.
  */
double Pucker_CP(const double* a1, const double* a2, const double* a3, 
                 const double* a4, const double* a5, const double* a6,
                 int N, double& amplitude, double& theta) 
{
  double dN = (double)N;
  double twopi_over_N = Constants::TWOPI / dN;
  Vec3 XYZ[6];
  XYZ[1].Assign( a1 );
  XYZ[2].Assign( a2 );
  XYZ[3].Assign( a3 );
  XYZ[4].Assign( a4 );
  if (N == 5) {
    XYZ[0].Assign( a5 ); // Ring apex
  } else if (N == 6) {
    XYZ[0].Assign( a6 ); // Ring apex
    XYZ[5].Assign( a5 );
  } else
    return -1.0; // Internal error
  // Calculate geometric center
  Vec3 Center(0.0, 0.0, 0.0);
  for (int i = 0; i < N; i++)
    Center += XYZ[i];
  Center /= dN;
  // Translate to center
  for (int i = 0; i < N; i++)
    XYZ[i] -= Center;
  // Calculate position vectors
  Vec3 R1(0.0, 0.0, 0.0);
  Vec3 R2(0.0, 0.0, 0.0);
  for (int i = 0; i < N; i++) {
    double factor = twopi_over_N * (double)i;
    double sin_val = sin( factor );
    double cos_val = cos( factor );
    R1[0] += XYZ[i][0] * sin_val;
    R2[0] += XYZ[i][0] * cos_val;
    R1[1] += XYZ[i][1] * sin_val;
    R2[1] += XYZ[i][1] * cos_val;
    R1[2] += XYZ[i][2] * sin_val; 
    R2[2] += XYZ[i][2] * cos_val;
  }
  // Calculate vector normal to plane
  Vec3 NXYZ = R1.Cross( R2 );
  // Normalize
  NXYZ.Normalize();
  // Calculate sums for determining q2/phi2
  // NOTE: Total amplitude Q should = sqrt(SUM[ Zn*Zn ])
  double sum1 = 0.0;
  double sum2 = 0.0;
  double Zn[6];
  for (int i = 0; i < N; i++) {
    double factor = 2.0 * twopi_over_N * (double)i;
    double cos_val = cos( factor );
    double sin_val = sin( factor ); 
    Zn[i] = XYZ[i] * NXYZ;
    sum1 += Zn[i] * cos_val;
    sum2 -= Zn[i] * sin_val;
  }
  // Calculate amplitude (q2, ==Q for N==5)
  double norm = sqrt(sum1*sum1 + sum2*sum2);
  amplitude = norm * sqrt( 2.0 / dN );
  // For even # coords (only 6 currently) calc extra pucker coord (q3)
  if (N == 6) {
    double q3 = 0.0;
    double mult = 1.0;
    for (int i = 0; i < N; i++) {
      q3 += mult * Zn[i]; // mult ~ pow( -1.0, i )
      mult = -mult;
    }
    q3 /= sqrt( dN );
    // Calculate theta
    theta = atan2( amplitude, q3 );
    // Update amplitude (Q for N==6)
    amplitude = sqrt( amplitude*amplitude + q3*q3 );
  }
  // Calculate pucker (phi2)
  double pucker = asin( sum2 / norm );
  if (sum1 < 0.0)
    pucker = Constants::PI - pucker;
  else if (pucker < 0.0)
    pucker += Constants::TWOPI;

  return pucker;
}

// CalcAngle()
double CalcAngle(const double* V1, const double* V2, const double* V3)
{
  double angle;
  double xij = V1[0] - V2[0];
  double yij = V1[1] - V2[1];
  double zij = V1[2] - V2[2];
  
  double xkj = V3[0] - V2[0];
  double ykj = V3[1] - V2[1];
  double zkj = V3[2] - V2[2];
  
  double rij = xij*xij + yij*yij + zij*zij;
  double rkj = xkj*xkj + ykj*ykj + zkj*zkj;

  if (rij > Constants::SMALL && rkj > Constants::SMALL) {
    angle = (xij*xkj + yij*ykj + zij*zkj) / sqrt(rij*rkj);
    if (angle > 1.0)
      angle = 1.0;
    else if (angle < -1.0)
      angle = -1.0;
    angle = acos(angle);
  } else
    angle = 0.0;

  return angle;
}
