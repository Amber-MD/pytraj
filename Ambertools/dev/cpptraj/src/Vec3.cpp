#include <cmath> // sqrt
#include "Vec3.h"
#include "CpptrajStdio.h"

/** Normalize vector. Return vector length. */
double Vec3::Normalize() {
  double r = sqrt( Magnitude2() );
  double b = 1.0 / r;
  V_[0] *= b; 
  V_[1] *= b; 
  V_[2] *= b;
  return r;
} 

double Vec3::Length() const { return sqrt( Magnitude2() ); }

void Vec3::Print(const char *Name) const {
  mprintf("    %s: %8.4f %8.4f %8.4f\n", Name, V_[0], V_[1], V_[2]);
}

/** Return the angle obtained from the dot product between this vector 
  * and U. Only works correctly if both are normalized beforehand.
  */
double Vec3::Angle(const Vec3& U) const {
  return acos( *this * U );
}

/** Return the angle obtained from the dot product between vectors V
  * and U, with sign determined from (VxU) dot Z. Assumes vectors 
  * are normalized.
  */
double Vec3::SignedAngle(const Vec3& U, const Vec3& Z) const {
  double dp = Angle( U );
  Vec3 Vec = Cross( U );
  double sign = Vec * Z;
  if (sign < 0) return -dp;
  return dp;
}
