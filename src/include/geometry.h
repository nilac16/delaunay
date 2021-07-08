#pragma once

#ifndef DELAUNAY_GEOMETRY_H
#define DELAUNAY_GEOMETRY_H


double dot_twovec(const double v1[static 2],
                  const double v2[static 2]);


/// \returns The magnitude of the vector resulting from the cross product of 
/// 2-vectors @p v1 and @p v2
///
double cross_twovec(const double v1[static 2],
                    const double v2[static 2]);


/// \param src
///     Vector to source point
/// \param dst
///     Vector to endpoint
/// \param v
///     The vector defined by @p dst - @p src
///
void relative_twovec(const double src[static 2],
                     const double dst[static 2],
                     double v[static 2]);


/// \returns | @p B - @p A |²
///
double segment_d_sqrd(const double A[static 2],
                      const double B[static 2]);


/// Finds the orientation of directed triangle ABC
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
