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

cpdef torsion(double[:, :, :] p):
   import numpy as np
   cdef double[:] out = np.empty(p.shape[0])
   cdef int i

   if p.shape[1] != 4 and p.shape[2] != 3:
       raise ValueError("shape of input array must be (n_frames, 4, 3)")

   for i in range(p.shape[0]):
       out[i] = math.degrees(C_Torsion(&p[i, 0, 0], &p[i, 1, 0],
                             &p[i, 2, 0], &p[i, 3, 0]))
   return np.asarray(out)

cpdef angle(double[:, :, :] p):
   import numpy as np
   cdef double[:] out = np.empty(p.shape[0])
   cdef int i

   if p.shape[1] != 3 and p.shape[2] != 3:
       raise ValueError("shape of input array must be (n_frames, 4, 3)")

   for i in range(p.shape[0]):
       out[i] = math.degrees(C_CalcAngle(&p[i, 0, 0], &p[i, 1, 0],
                             &p[i, 2, 0]))
   return np.asarray(out)
