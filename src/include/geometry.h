#pragma once

#ifndef DELAUNAY_GEOMETRY_H
#define DELAUNAY_GEOMETRY_H

/** TODO:
 *  Go through this file and remove forward decls for unused functions
 *  Make them static if they are used internally, otherwise, erase the 
 *  definitions as well
 * 
 *  ALSO:
 *  Watch for extreme UB from some of these restricted pointers, especially
 *  the function collinear_center()
**/

double dot_twovec(const double v1[static 2],
                  const double v2[static 2]);

double cross_twovec(const double v1[static 2],
                    const double v2[static 2]);

void relative_twovec(const double src[static 2],
                     const double dst[static 2],
                     double v[static 2]);


/// Returns | B - A |²
///
double segment_d_sqrd(const double A[static 2],
                      const double B[static 2]);


/// Returns the orientation of directed triangle ABC
///
/// \returns 1 if the triangle is positively oriented, -1 if negatively 
///     oriented, or 0 if degenerate (A, B, and C are collinear)
///
int triangle_orientation(const double A[static restrict 2],
                         const double B[static restrict 2],
                         const double C[static restrict 2]);


/// Finds the center of the circumcircle of triangle ABC, and stores it in
/// @p r. The third value is the radius squared
///
void triangle_circumcircle(const double A[static restrict 2],
                           const double B[static restrict 2],
                           const double C[static restrict 2],
                           double r[static restrict 3]);


/// Finds the index of the node at the center of the three passed collinear
/// nodes
///
int collinear_center(const double *restrict rs[static 3]);


#endif //DELAUNAY_GEOMETRY_H
