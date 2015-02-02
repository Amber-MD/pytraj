#ifndef INC_DISTROUTINES_H
#define INC_DISTROUTINES_H
#include "Box.h"
/*! \file DistRoutines.h
    \brief A collection of routines used to calculate distance.
 */
// TODO: Move ImagingType to Box, include Box instead of Vec3, create recipInfo type
/// Potential imaging types 
enum ImagingType { NOIMAGE=0, ORTHO, NONORTHO };
double DIST2_ImageNonOrtho(Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&);
Vec3 MinImagedVec(Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&);
double DIST2_ImageNonOrthoRecip(Vec3 const&, Vec3 const&, double, int*, Matrix_3x3 const&);
double DIST2_ImageOrtho(Vec3 const&, Vec3 const&, Box const&);
double DIST2_NoImage(const double*, const double*);
double DIST2_NoImage( Vec3 const&, Vec3 const& );
double DIST_NoImage( Vec3 const&, Vec3 const& );
double DIST2(const double*, const double*, ImagingType, Box const&, 
             Matrix_3x3 const&, Matrix_3x3 const&);
#endif
