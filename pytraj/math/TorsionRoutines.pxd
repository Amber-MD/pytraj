# distutil: language = c++

cdef extern from "TorsionRoutines.h":
    double Torsion(const double *, const double *, const double *, const double *)
    double Pucker_AS(const double*, const double*, const double*, const double*, 
                             const double*, double&)
    double Pucker_CP(const double*, const double*, const double*, const double*, 
                             const double*, const double*, int, double&, double&)
    double CalcAngle(const double*, const double*, const double*)
