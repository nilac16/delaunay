/**
 *  Different geometric primitives needed for triangulation. This probably 
 *  doesn't need to be its own file, but the file that writes the triangulation
 *  to an image on disk uses non-static functions exposed in the header, so I 
 *  will not condense this into delaunay.c
 * 
 *  This uses x86 vector intrinsics by default. To disable this behavior, 
 *  define NO_INTRINSICS. The simplest way to do so is probably just to compile
 *  with -DNO_INTRINSICS
 */
#include <math.h>
#include <stdio.h>
#include <x86intrin.h>
#include "include/geometry.h"

#ifdef __GNUC__
#define gnu_attribute(...) __attribute__((__VA_ARGS__))
#else
#define gnu_attribute(...)
#endif //__GNUC__


gnu_attribute(nonnull, pure)
double dot_twovec(const double *restrict v1,
                  const double *restrict v2)
{
    #ifdef NO_INTRINSICS
    return v1[0] * v2[0] + v1[1] * v2[1];
    #else
    __m128d vec1 = _mm_load_pd(v1);
    __m128d vec2 = _mm_load_pd(v2);
    vec1 = _mm_dp_pd(vec1, vec2, 0x31);
    double ret;
    _mm_store_sd(&ret, vec1);
    return ret;
    #endif //NO_INTRINSICS
}

gnu_attribute(const)
/// Crossing two 2D vectors can be likened to a dot product, where one vector 
/// is transposed (shuffled) and its lower (upper bits!) are negated (the 
/// blend with out_n). This returns a scalar, the magnitude of the normal 
/// vector to the plane defined by @p v1 and @p v2
///
static double cross_twovec_xmm(const __m128d v1, const __m128d v2)
{
    __m128d out = _mm_shuffle_pd(v2, v2, _MM_SHUFFLE2(0, 1));
    __m128d out_n = _mm_sub_pd(_mm_setzero_pd(), out);
    out = _mm_blend_pd(out, out_n, 0x2);
    out = _mm_dp_pd(v1, out, 0x31);
    double ret;
    _mm_store_sd(&ret, out);
    return ret;
}

gnu_attribute(nonnull, pure)
double cross_twovec(const double *restrict v1,
                    const double *restrict v2)
{
    #ifdef NO_INTRINSICS
    return v1[0] * v2[1] - v1[1] * v2[0];
    #else
    __m128d vec1 = _mm_load_pd(v1);
    __m128d vec2 = _mm_load_pd(v2);
    return cross_twovec_xmm(vec1, vec2);
    #endif //NO_INTRINSICS
}

gnu_attribute(nonnull)
/// Stores the 2-vector from @p src to @p dst in @p v
///
/// \param src
///     Vector to start point
/// \param dst
///     Vector endpoint
/// \param v
///     Buffer to write the vector to
///
void relative_twovec(const double *restrict src,
                     const double *restrict dst,
                     double *restrict v)
{
    #ifdef NO_INTRINSICS
    v[0] = dst[0] - src[0];
    v[1] = dst[1] - src[1];
    #else
    __m128d srcv = _mm_load_pd(src);
    __m128d dstv = _mm_load_pd(dst);
    srcv = _mm_sub_pd(dstv, srcv);
    _mm_store_pd(v, srcv);
    #endif //NO_INTRINSICS
}

gnu_attribute(nonnull, pure)
double segment_d_sqrd(const double *restrict A, const double *restrict B)
{
    #ifdef NO_INTRINSICS
    double dr[2];
    relative_twovec(A, B, dr);
    return dr[0] * dr[0] + dr[1] * dr[1];
    #else
    __m128d vA = _mm_load_pd(A);
    __m128d vB = _mm_load_pd(B);
    __m128d dv = _mm_sub_pd(vB, vA);
    dv = _mm_dp_pd(dv, dv, 0x31);
    double out;
    _mm_store_sd(&out, dv);
    return out;
    #endif //NO_INTRINSICS
}

gnu_attribute(const)
static int signum(double x)
{
    if (x < 0) {
        return -1;
    } else {
        return x > 0;
    }
}

gnu_attribute(nonnull, pure, hot)
int triangle_orientation(const double *restrict A,
                         const double *restrict B,
                         const double *restrict C)
{
    double orientation;
    #ifdef NO_INTRINSICS
    double AB[2];
    double BC[2];
    relative_twovec(A, B, AB);
    relative_twovec(B, C, BC);
    orientation = cross_twovec(AB, BC);
    #else
    __m128d AB = _mm_sub_pd(_mm_load_pd(B), _mm_load_pd(A));
    __m128d BC = _mm_sub_pd(_mm_load_pd(C), _mm_load_pd(B));
    orientation = cross_twovec_xmm(AB, BC);
    #endif //NO_INTRINSICS
    return signum(orientation);
}

gnu_attribute(nonnull)
void triangle_circumcircle(const double *restrict A,
                           const double *restrict B,
                           const double *restrict C,
                           double *restrict r)
{
    //            1      | (B² - C²)(By - Ay) + (B² - A²)(Cy - By) |
    //  r =  ----------- |                                         |
    //       2 |AB × BC| | (C² - B²)(Bx - Ax) + (A² - B²)(Cx - Bx) |
    //
    #ifdef NO_INTRINSICS
    double AB[2], BC[2];
    relative_twovec(A, B, AB);
    relative_twovec(B, C, BC);
    double denom = 2 * cross_twovec(AB, BC);
    double B_sqrd = dot_twovec(B, B);
    double terms[2] = {
        dot_twovec(C, C) - B_sqrd,
        dot_twovec(A, A) - B_sqrd
    };
    r[0] = -(terms[0] * AB[1] + terms[1] * BC[1]) / denom;
    r[1] = (terms[0] * AB[0] + terms[1] * BC[0]) / denom;
    #else
    __m128d v1 = _mm_load_pd(B);
    __m128d AB = _mm_sub_pd(v1, _mm_load_pd(A));
    __m128d BC = _mm_sub_pd(_mm_load_pd(C), v1);
    __m128d denom = _mm_mul_pd(_mm_set1_pd(2.0), _mm_set1_pd(cross_twovec_xmm(AB, BC)));
    // AB = < Bx - Ax, By - Ay >
    // AB = < Cx - Bx, Cy - By >
    // v1 = < B, B >
    // denom = < |AB × BC|, |AB × BC| >
    v1 = _mm_dp_pd(v1, v1, 0x33);
    // v1 = < B², B² >
    __m128d v2 = _mm_load_pd(C);
    // v2 = < Cx, Cy >
    AB = _mm_mul_pd(_mm_sub_pd(_mm_dp_pd(v2, v2, 0x33), v1), AB);
    // AB = < (C² - B²)ABx, (C² - B²)ABy >   C²        - B²
    v2 = _mm_load_pd(A);
    // v2 = < Ax, Ay >
    BC = _mm_mul_pd(_mm_sub_pd(_mm_dp_pd(v2, v2, 0x33), v1), BC);
    // BC = < (A² - B²)BCx, (A² - B²)BCy >
    v1 = _mm_add_pd(AB, BC);
    // v1 = < (C² - B²)ABx + (A² - B²)BCx, (C² - B²)ABy + (A² - B²)BCy >
    v1 = _mm_div_pd(v1, denom);
    v1 = _mm_shuffle_pd(v1, v1, _MM_SHUFFLE2(0, 1));
    v2 = _mm_sub_pd(_mm_setzero_pd(), v1);
    // negate only the first component of v1
    v1 = _mm_blend_pd(v1, v2, 0x1);
    _mm_store_pd(r, v1);
    #endif //NO_INTRINSICS
    r[2] = segment_d_sqrd(r, A);
}

/* gnu_attribute(nonnull)
int in_circumcircle(const double *restrict A, const double *restrict B,
                    const double *restrict C, const double *restrict D)
{
    double A_comp = A[1];
    double minor0, minor1, minor2;
    {
        const double A_sqrd = dot_twovec(A, A);
        const double B_sqrd = dot_twovec(B, B);
        const double C_sqrd = dot_twovec(C, C);
        const double D_sqrd = dot_twovec(D, D);
        minor0 = (B[1] - A_comp) * (C_sqrd - A_sqrd)
               - (C[1] - A_comp) * (B_sqrd - A_sqrd);
        minor1 = (C[1] - A_comp) * (D_sqrd - A_sqrd)
               - (D[1] - A_comp) * (C_sqrd - A_sqrd);
        minor2 = (D[1] - A_comp) * (B_sqrd - A_sqrd)
               - (B[1] - A_comp) * (D_sqrd - A_sqrd);
    }
    A_comp = A[0];
    const double det = (D[0] - A_comp) * minor0
                     + (B[0] - A_comp) * minor1
                     + (C[0] - A_comp) * minor2;
    return det < 0;
} */

gnu_attribute(nonnull)
int collinear_center(const double *restrict rs[static 3])
{
    double AB[2], AC[2];
    relative_twovec(rs[2], rs[0], AB);
    relative_twovec(rs[2], rs[1], AC);
    if (dot_twovec(AB, AC) < 0) {
        return 2;
    } else {
        relative_twovec(rs[1], rs[0], AB);
        relative_twovec(rs[1], rs[2], AC);
        return dot_twovec(AB, AC) < 0;
    }
}
