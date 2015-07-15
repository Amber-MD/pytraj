# distutil: language = c++
import math

cdef extern from "TorsionRoutines.h" nogil:
    # create alias to avoid: ambiguous overloaded method
    double C_Torsion "Torsion" (const double *, const double *, const double *, const double *)
    double Pucker_AS(const double*, const double*, const double*, const double*, 
                             const double*, double&)
    double Pucker_CP(const double*, const double*, const double*, const double*, 
                             const double*, const double*, int, double&, double&)
    double C_CalcAngle "CalcAngle" (const double*, const double*, const double*)

cpdef torsion(double[:, :] p):
   return math.degrees(C_Torsion(&p[0, 0], &p[1, 0], &p[2, 0], &p[3, 0]))

cpdef angle(double[:, :] p):
   return math.degrees(C_CalcAngle(&p[0, 0], &p[1, 0], &p[2, 0]))
