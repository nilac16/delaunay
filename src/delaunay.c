#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "include/geometry.h"
#include "include/delaunay.h"

#ifndef __WIN32
#include <pthread.h>
#else
#include <windows.h>
#endif

#ifdef __GNUC__
#define gnu_attribute(...) __attribute__((__VA_ARGS__))
#else
#define gnu_attribute(...)
#endif //__GNUC__





#define DIM_TOGGLE(dim) (!dim)

enum {
    DIM_X = 0,
    DIM_Y = 1
};


/// Expands to the necessary calculation for a counterclockwise rotation,
/// i times, from starting orientation o. This assumes that i will never be 
/// less than -3 (all calls within this program have i either 1 or -1, see the 
/// enum below this)
///
#define ROT(o, i) ((3 + o + i) % 3)

/// In order to go clockwise around a point, we must rotate each triangle to 
/// its counterclockwise orientation before moving to the next, and vice versa
///
enum {
    ROT_CW  = -1,
    ROT_CCW =  1
};

#define DIR_TOGGLE(dir) (-dir)


/// The orientation (rot) of a manifold triangle that points in the labeled 
/// direction. This is only enforced for manifold triangles, so that traversing 
/// can be fast, and not require maintaining a separate integer for their 
/// orientations
///
enum {
    MANIFOLD_CW  = 0,   // Clockwise around convex hull
    MANIFOLD_CCW = 1,   // Counterclockwise around convex hull
    MANIFOLD_IN  = 2    // Inwards, into the triangulation
};

#define NEXT_MANIFOLD(tri, dir) (tri)->adjacents[dir].next


#define TRI_VERTEX(tri, rot) (tri)->vertices[rot]


gnu_attribute(malloc)
static struct delaunay_triangle_pool *new_triangle_pool(size_t N)
{
    struct delaunay_triangle_pool *p = malloc(sizeof *p + (N * (sizeof *p->start)));
    if (!p) {
        return NULL;
    }
    p->end = p->start;
    return p;
}

gnu_attribute(nonnull)
void free_triangle_pool(struct delaunay_triangle_pool *p)
{
    free(p);
}


#ifdef __WIN32
static CRITICAL_SECTION crit;
#endif //__WIN32


gnu_attribute(nonnull, malloc)
/// Reserves @p N triangles on the stack allocator at @p p. This is the only 
/// critical section in the entire algorithm, so all of the concurrency 
/// primitives are clustered here
///
static struct delaunay_triangle *reserve_tris(struct delaunay_triangle_pool *p,
                                              size_t N)
{
    #ifndef __WIN32
    static pthread_mutex_t mutex;
    pthread_mutex_lock(&mutex);
    #else
    EnterCriticalSection(&crit);
    #endif //__WIN32

    struct delaunay_triangle *t = p->end;
    p->end += N;

    #ifndef __WIN32
    pthread_mutex_unlock(&mutex);
    #else
    LeaveCriticalSection(&crit);
    #endif //__WIN32

    return t;
}

gnu_attribute(nonnull, hot)
/// Updates the corresponding triangles' sides to point at each other
///
static void retarget_tri(struct delaunay_triangle *restrict t1, int o1,
                         struct delaunay_triangle *restrict t2, int o2)
{
    t1->adjacents[o1].next = t2;
    t1->adjacents[o1].next_side = o2;
    t2->adjacents[o2].next = t1;
    t2->adjacents[o2].next_side = o1;
}

gnu_attribute(nonnull, malloc)
/// A line segment using Shewchuk's triangle structure is composed of two 
/// manifold triangles. Here, I make the manifold node index 2, and it is 
/// shared by each triangle. Both triangles have the same orientation
///
static struct delaunay_triangle *make_segment(struct delaunay_triangle_pool *p,
                                              const double *restrict v0,
                                              const double *restrict v1)
{
    struct delaunay_triangle *ts = reserve_tris(p, 2);
    TRI_VERTEX(ts, MANIFOLD_CW) = v0;
    TRI_VERTEX(ts, MANIFOLD_CCW) = v1;
    TRI_VERTEX(ts, MANIFOLD_IN) = NULL;

    TRI_VERTEX(ts + 1, MANIFOLD_CW) = v1;
    TRI_VERTEX(ts + 1, MANIFOLD_CCW) = v0;
    TRI_VERTEX(ts + 1, MANIFOLD_IN) = NULL;

    retarget_tri(ts, MANIFOLD_IN, ts + 1, MANIFOLD_IN);
    retarget_tri(ts, MANIFOLD_CW, ts + 1, MANIFOLD_CCW);
    retarget_tri(ts + 1, MANIFOLD_CW, ts, MANIFOLD_CCW);

    return ts;
}

gnu_attribute(nonnull)
/// Checks the linearity of the triangle formed by @p nodes. If the nodes are
/// in clockwise order, flips them to CCW. If they are collinear, returns a 
/// sentinel value
///
static int sort_ccw(const double *restrict v0,
                    const double *restrict *v1,
                    const double *restrict *v2)
{
    int o = triangle_orientation(v0, *v1, *v2);
    if (o == -1) {
        const double *tmp = *v2;
        *v2 = *v1;
        *v1 = tmp;
    }
    return !o;
}

gnu_attribute(nonnull, cold)
/// In the case of collinear points, it is not possible to make a real 
/// triangle, so we must make two attached line segments instead
///
static struct delaunay_triangle *double_segments(struct delaunay_triangle *ts,
                                                 const double *vs[static 3])
{
    //puts(" DEGENERATE TRIANGLE");
    {
        int center = collinear_center(vs);
        if (center != 1) {
            const double *tmp = vs[center];
            vs[center] = vs[1];
            vs[center] = tmp;
        }
    }
    TRI_VERTEX(ts, MANIFOLD_CW) = vs[0];
    TRI_VERTEX(ts, MANIFOLD_CCW) = vs[1];
    TRI_VERTEX(ts, MANIFOLD_IN) = NULL;

    TRI_VERTEX(ts + 1, MANIFOLD_CW) = vs[1];
    TRI_VERTEX(ts + 1, MANIFOLD_CCW) = vs[0];
    TRI_VERTEX(ts + 1, MANIFOLD_IN) = NULL;

    TRI_VERTEX(ts + 2, MANIFOLD_CW) = vs[2];
    TRI_VERTEX(ts + 2, MANIFOLD_CCW) = vs[1];
    TRI_VERTEX(ts + 2, MANIFOLD_IN) = NULL;

    TRI_VERTEX(ts + 3, MANIFOLD_CW) = vs[1];
    TRI_VERTEX(ts + 3, MANIFOLD_CCW) = vs[2];
    TRI_VERTEX(ts + 3, MANIFOLD_IN) = NULL;

    retarget_tri(ts, MANIFOLD_IN, ts + 1, MANIFOLD_IN);
    retarget_tri(ts + 2, MANIFOLD_IN, ts + 3, MANIFOLD_IN);
    
    retarget_tri(ts, MANIFOLD_CCW, ts + 1, MANIFOLD_CW);
    retarget_tri(ts + 1, MANIFOLD_CCW, ts + 2, MANIFOLD_CW);
    retarget_tri(ts + 2, MANIFOLD_CCW, ts + 3, MANIFOLD_CW);
    retarget_tri(ts + 3, MANIFOLD_CCW, ts, MANIFOLD_CW);

    return ts;
}

gnu_attribute(nonnull, malloc)
static struct delaunay_triangle *make_triangle(struct delaunay_triangle_pool *p,
                                               const double *restrict v0,
                                               const double *restrict v1,
                                               const double *restrict v2)
{
    struct delaunay_triangle *ts = reserve_tris(p, 4);
    if (sort_ccw(v0, &v1, &v2)) {
        const double *vs[3] = { v0, v1, v2 };
        return double_segments(ts, vs);
    }
    /// TODO: There is probably a more efficient way to do this
    TRI_VERTEX(ts, 0) = v0;
    TRI_VERTEX(ts, 1) = v1;
    TRI_VERTEX(ts, 2) = v2;

    TRI_VERTEX(ts + 1, MANIFOLD_CW) = v2;
    TRI_VERTEX(ts + 1, MANIFOLD_CCW) = v1;
    TRI_VERTEX(ts + 1, MANIFOLD_IN) = NULL;

    TRI_VERTEX(ts + 2, MANIFOLD_CW) = v0;
    TRI_VERTEX(ts + 2, MANIFOLD_CCW) = v2;
    TRI_VERTEX(ts + 2, MANIFOLD_IN) = NULL;

    TRI_VERTEX(ts + 3, MANIFOLD_CW) = v1;
    TRI_VERTEX(ts + 3, MANIFOLD_CCW) = v0;
    TRI_VERTEX(ts + 3, MANIFOLD_IN) = NULL;

    retarget_tri(ts, 0, ts + 1, MANIFOLD_IN);
    retarget_tri(ts, 1, ts + 2, MANIFOLD_IN);
    retarget_tri(ts, 2, ts + 3, MANIFOLD_IN);

    retarget_tri(ts + 1, MANIFOLD_CCW, ts + 2, MANIFOLD_CW);
    retarget_tri(ts + 2, MANIFOLD_CCW, ts + 3, MANIFOLD_CW);
    retarget_tri(ts + 3, MANIFOLD_CCW, ts + 1, MANIFOLD_CW);

    triangle_circumcircle(v0, v1, v2, ts->circumcenter);

    return ts + 1;
}

gnu_attribute(nonnull, hot)
/// Moves the triangle pointer @p t to the triangle neighboring it at @p rot.
/// @p rot is then updated to the neighboring triangle's orientation, such 
/// that successive calls to this would oscillate between the two triangles
///
static void next_triangle(struct delaunay_triangle *restrict *restrict t, int *rot)
{
    int rot_new = (*t)->adjacents[*rot].next_side;
    *t = (*t)->adjacents[*rot].next;
    *rot = rot_new;
}

gnu_attribute(nonnull, hot)
static void next_rotate(struct delaunay_triangle *restrict *restrict t, int *rot, int direction)
{
    next_triangle(t, rot);
    *rot = ROT(*rot, direction);
}

gnu_attribute(nonnull)
/// Flips the edge associated with triangle @p t at orientation @p rot 
/// The triangle that is passed to this retains its vertex that is NOT in 
/// @p direction! Use this when deleting edges while merging!
///
static void flip_edge(struct delaunay_triangle *t, int rot, int direction)
{
    struct delaunay_triangle *t2 = t;
    int rot2 = rot;
    next_triangle(&t2, &rot2);
    TRI_VERTEX(t, ROT(rot, direction)) = TRI_VERTEX(t2, rot2);
    TRI_VERTEX(t2, ROT(rot2, direction)) = TRI_VERTEX(t, rot);
    direction = DIR_TOGGLE(direction);
    retarget_tri(t, rot, t2->adjacents[ROT(rot2, direction)].next, t2->adjacents[ROT(rot2, direction)].next_side);
    retarget_tri(t2, rot2, t->adjacents[ROT(rot, direction)].next, t->adjacents[ROT(rot, direction)].next_side);
    retarget_tri(t, ROT(rot, direction), t2, ROT(rot2, direction));
}

gnu_attribute(nonnull, pure)
/// Searches around the convex hull of @p t (which must be a manifold triangle) 
/// until it finds an inflection point with positive second derivative. 
/// Scanning occurs in @p direction around the hull; use MANIFOLD_CCW for the 
/// left triangulation and MANIFOLD_CW for the right triangulation. This 
/// ensures that the result triangles are ABOVE the nascent base LR edge
///
static struct delaunay_triangle *minimum_node(struct delaunay_triangle *t,
                                              int direction, int dim)
{
    // scan until it inflects down
    while (TRI_VERTEX(t, direction)[dim] < TRI_VERTEX(t, !direction)[dim]) {
        t = NEXT_MANIFOLD(t, direction);
    }
    // scan until it inflects up
    while (TRI_VERTEX(t, direction)[dim] > TRI_VERTEX(t, !direction)[dim]) {
        t = NEXT_MANIFOLD(t, direction);
    }
    return t;
}

gnu_attribute(nonnull)
/// Determines if the (directed) triangles defined by @p L -> @p R and both of
/// the neighbors of ( @p t, @p rot ) are both positively oriented. If this is 
/// the case, then the entire convex hull that @p t is on is to the LEFT of 
/// directed edge L->R
///
static int lower_tangent(const double *restrict L, const double *restrict R,
                         struct delaunay_triangle *t, int rot)
{
    int o = triangle_orientation(L, R, TRI_VERTEX(t, !rot));
    if (o < 0) {
        return 0;
    }
    o = triangle_orientation(L, R, TRI_VERTEX(NEXT_MANIFOLD(t, !rot), rot));
    return !(o < 0);
}

gnu_attribute(nonnull)
/// Computes the lower common tangent of both convex hulls at @p tleft and 
/// @p tright. 
/// TODO: Reimplement A&W's common tangent algorithm
///
/// \param tleft
///     Reference to left triangle pointer. On entry, this must be a triangle 
///     on the L triangulation's manifold. On exit, this is the triangle 
///     containing the L vertex of the base LR edge, at its counterclockwise 
///     index
/// \param tright
///     Reference to right triangle pointer. The same entry restrictions 
///     apply. On exit, this contains the R vertex at its clockwise index
/// \param dim
///     The current kd-tree cutting dimension. This is needed to minimize the 
///     initial nodes passed to A&W's algorithm
///
static void find_base_LR_edge(struct delaunay_triangle *restrict *restrict tleft,
                              struct delaunay_triangle *restrict *restrict tright,
                              int dim)
{
    *tleft = minimum_node(*tleft, MANIFOLD_CCW, !dim);
    *tright = minimum_node(*tright, MANIFOLD_CW, !dim);
    const double *restrict L = TRI_VERTEX(*tleft, MANIFOLD_CCW);
    const double *restrict R = TRI_VERTEX(*tright, MANIFOLD_CW);
    // iterate on the right until finding the lower tangent
    // then make sure that the left point is the tangent
    // if not, move it and reiterate
    while (1) {
        while (!lower_tangent(L, R, *tright, MANIFOLD_CW)) {
            *tright = NEXT_MANIFOLD(*tright, MANIFOLD_CCW);
            R = TRI_VERTEX(*tright, MANIFOLD_CW);
        }
        if (!lower_tangent(L, R, *tleft, MANIFOLD_CCW)) {
            *tleft = NEXT_MANIFOLD(*tleft, MANIFOLD_CW);
            L = TRI_VERTEX(*tleft, MANIFOLD_CCW);
        } else {
            break;
        }
    }
}

gnu_attribute(nonnull, malloc)
/// Connects the manifold triangles @p t1 and @p t2 at their vertices 
/// specified by  @p o1 and @p o2, respectively. The new edge is a closed 
/// segment that is spliced into the triangulations at @p t1 and @p t2. 
///
static struct delaunay_triangle *
make_edge(struct delaunay_triangle_pool *p,
          struct delaunay_triangle *restrict t1, int o1,
          struct delaunay_triangle *restrict t2, int o2)
{
    struct delaunay_triangle *ts = reserve_tris(p, 2);

    TRI_VERTEX(ts, MANIFOLD_CW) = TRI_VERTEX(t2, o2);
    TRI_VERTEX(ts, MANIFOLD_CCW) = TRI_VERTEX(t1, o1);
    TRI_VERTEX(ts, MANIFOLD_IN) = NULL;

    TRI_VERTEX(ts + 1, MANIFOLD_CW) = TRI_VERTEX(t1, o1);
    TRI_VERTEX(ts + 1, MANIFOLD_CCW) = TRI_VERTEX(t2, o2);
    TRI_VERTEX(ts + 1, MANIFOLD_IN) = NULL;

    retarget_tri(ts, MANIFOLD_IN, ts + 1, MANIFOLD_IN);

    retarget_tri(ts, MANIFOLD_CW, t1->adjacents[!o1].next, t1->adjacents[!o1].next_side);
    retarget_tri(ts, MANIFOLD_CCW, t2->adjacents[!o2].next, t2->adjacents[!o2].next_side);

    retarget_tri(ts + 1, MANIFOLD_CW, t2, !o2);
    retarget_tri(ts + 1, MANIFOLD_CCW, t1, !o1);

    return ts + 1;
}

gnu_attribute(nonnull, pure, hot)
/// Checks if the passed edge triangle's candidate note intersects with the 
/// next candidate around it
///
/// \param cand
///     Manifold triangle bordering the start tri
/// \param direction
///     Direction to rotate around the LR point
///
static int candidate_collision(struct delaunay_triangle *cand,
                               double circum[static 3], int direction)
{
    direction = DIR_TOGGLE(direction);
    int rot = MANIFOLD_IN;
    next_rotate(&cand, &rot, direction);
    // but we have to go again, because this node only defines the original
    // candidate node
    next_rotate(&cand, &rot, direction);
    return TRI_VERTEX(cand, rot) && (segment_d_sqrd(TRI_VERTEX(cand, rot), circum) < circum[2]);
}

gnu_attribute(nonnull)
static void swap_sides(struct delaunay_triangle *t, int o1, int o2)
{
    const double *v = TRI_VERTEX(t, o1);
    struct delaunay_triangle *n = t->adjacents[o1].next;
    int x = t->adjacents[o1].next_side;
    TRI_VERTEX(t, o1) = TRI_VERTEX(t, o2);
    retarget_tri(t, o1, t->adjacents[o2].next, t->adjacents[o2].next_side);
    TRI_VERTEX(t, o2) = v;
    retarget_tri(t, o2, n, x);
}

gnu_attribute(nonnull)
/// Ensures that the manifold pointer of triangle @p t is at index 2. If it 
/// is not, multiple pointer shuffles occur to correct this
///
/// \param t
///     Manifold triangle to validate
///
static void validate_manifold(struct delaunay_triangle *t)
{
    int i = 0;
    while (TRI_VERTEX(t, i)) {
        i++;
    }
    if (i == MANIFOLD_IN) {
        return;
    }
    swap_sides(t, i, MANIFOLD_IN);
    swap_sides(t, i, !i);
}

gnu_attribute(nonnull)
/// Finds a candidate that doesn't clash with the next candidate in the list. 
/// Goes @p direction around the list
///
/// \param tedge
///     This should be either tleft or tright
/// \param rot
///     The location of the candidate node on @p tedge
/// \param direction
///     The direction to rotate around the LR vertex
///
static int submit_candidate(struct delaunay_triangle *tedge, int rot,
                            int direction,
                            const double *restrict vl,
                            const double *restrict vr)
{
    double circum[3];
    while (1) {
        if (triangle_orientation(vl, vr, TRI_VERTEX(tedge, rot)) < 0) {
            return 0;
        }
        triangle_circumcircle(vl, vr, TRI_VERTEX(tedge, rot), circum);
        if (candidate_collision(tedge, circum, direction)) {
            // flip tedge and try again
            struct delaunay_triangle *tother = NEXT_MANIFOLD(tedge, MANIFOLD_IN);
            flip_edge(tedge, MANIFOLD_IN, direction);
            validate_manifold(tother);
        } else {
            break;
        }
    }
    return 1;
}

gnu_attribute(nonnull)
static struct delaunay_triangle *
merge_triangulations(struct delaunay_triangle_pool *p,
                     struct delaunay_triangle *restrict tleft,
                     struct delaunay_triangle *restrict tright,
                     int dim)
{
    find_base_LR_edge(&tleft, &tright, dim);
    const double *restrict vl, *restrict vr;
    struct delaunay_triangle *ts = make_edge(p, tleft, MANIFOLD_CCW, tright, MANIFOLD_CW);
    // ts now points to the upper manifold triangle of the base LR edge
merge_find_candidates:
    vl = TRI_VERTEX(tleft, MANIFOLD_CCW);
    vr = TRI_VERTEX(tright, MANIFOLD_CW);
    int l = submit_candidate(tleft, MANIFOLD_CW, ROT_CCW, vl, vr);
    int r = submit_candidate(tright, MANIFOLD_CCW, ROT_CW, vl, vr);
    if (l) {
        if (r) {
            double circum[3];
            triangle_circumcircle(vl, vr, TRI_VERTEX(tleft, MANIFOLD_CW), circum);
            if (segment_d_sqrd(circum, TRI_VERTEX(tright, MANIFOLD_CCW)) < circum[2]) {
                goto merge_accept_right;
            } else {
                goto merge_accept_left;
            }
        } else {
            goto merge_accept_left;
        }
    } else {
        if (r) {
            goto merge_accept_right;
        } else {
            return ts;
        }
    }
merge_accept_right:
    flip_edge(ts, MANIFOLD_CW, ROT_CCW);
    triangle_circumcircle(TRI_VERTEX(tright, 0), TRI_VERTEX(tright, 1), TRI_VERTEX(tright, 2), tright->circumcenter);
    tright = NEXT_MANIFOLD(ts, MANIFOLD_CW);
    goto merge_find_candidates;
merge_accept_left:
    flip_edge(ts, MANIFOLD_CCW, ROT_CW);
    triangle_circumcircle(TRI_VERTEX(tleft, 0), TRI_VERTEX(tleft, 1), TRI_VERTEX(tleft, 2), tleft->circumcenter);
    tleft = NEXT_MANIFOLD(ts, MANIFOLD_CCW);
    goto merge_find_candidates;
}

gnu_attribute(nonnull, pure, hot)
static int xcmp(const void *restrict v1, const void *restrict v2)
{
    double x1 = **(const double **)v1;
    double x2 = **(const double **)v2;
    return (x1 > x2) - (x1 < x2);
}

gnu_attribute(nonnull, pure, hot)
static int ycmp(const void *restrict v1, const void *restrict v2)
{
    double y1 = *(*(const double **)v1 + 1);
    double y2 = *(*(const double **)v2 + 1);
    return (y1 < y2) - (y1 > y2);
}

gnu_attribute(nonnull)
/// Serial version
///
static struct delaunay_triangle *build_2dtree_s(struct delaunay_triangle_pool *p,
                                                size_t N,
                                                const double *refs[static N],
                                                int dim)
{
    static int (*const cmps[])(const void *, const void *) = {
        xcmp,
        ycmp
    };
    if (N <= 3) {
        if (N == 3) {
            return make_triangle(p, refs[0], refs[1], refs[2]);
        } else {
            return make_segment(p, refs[0], refs[1]);
        }
    }
    struct delaunay_triangle *restrict t_left, *restrict t_right;
    qsort(refs, N, sizeof *refs, cmps[dim]);
    size_t med = N / 2;
    t_left = build_2dtree_s(p, med, refs, DIM_TOGGLE(dim));
    t_right = build_2dtree_s(p, (N + 1) / 2, refs + med, DIM_TOGGLE(dim));
    t_left = merge_triangulations(p, t_left, t_right, dim);
    return t_left;
}


#ifndef __WIN32
static void *delaunay_entry(void *args);
#else
DWORD delaunay_entry(void *args);
#endif //__WIN32


struct delaunay_args {
    struct delaunay_triangle_pool *p;
    size_t N;
    const double **refs;
    int dim;
    unsigned int n_threads;

    #ifdef __WIN32
    struct delaunay_triangle *result;
    #endif //__WIN32
};

gnu_attribute(nonnull)
/// Implicitly builds a 2d-tree, then tears it down
///
static struct delaunay_triangle *build_2dtree(struct delaunay_triangle_pool *p,
                                              size_t N,
                                              const double *refs[static N],
                                              int dim, unsigned int n_threads)
{
    static int (*const cmps[])(const void *, const void *) = {
        xcmp,
        ycmp
    };
    if (n_threads == 1) {
        return build_2dtree_s(p, N, refs, dim);
    }
    if (N <= 3) {
        if (N == 3) {
            return make_triangle(p, refs[0], refs[1], refs[2]);
        } else {
            return make_segment(p, refs[0], refs[1]);
        }
    }
    struct delaunay_triangle *restrict t_left, *restrict t_right;
    qsort(refs, N, sizeof *refs, cmps[dim]);
    size_t med = N / 2;

    struct delaunay_args args;
    args.p = p;
    args.N = (N + 1) / 2;
    args.refs = refs + med;
    args.dim = DIM_TOGGLE(dim);
    args.n_threads = n_threads / 2;

    #ifndef __WIN32
    pthread_t r;
    pthread_create(&r, NULL, delaunay_entry, &args);
    #else
    HANDLE r;
    DWORD hThread;
    r = CreateThread(NULL, 0, delaunay_entry, &args, 0, &hThread);
    #endif //__WIN32

    t_left = build_2dtree(p, med, refs, DIM_TOGGLE(dim), (n_threads + 1) / 2);

    #ifndef __WIN32
    pthread_join(r, (void **)&t_right);
    #else
    WaitForSingleObject(r, INFINITE);
    t_right = args.result;
    #endif //__WIN32

    t_left = merge_triangulations(p, t_left, t_right, dim);
    return t_left;
}

#ifndef __WIN32
gnu_attribute(nonnull, malloc)
static void *delaunay_entry(void *args)
{
    return build_2dtree(((struct delaunay_args *)args)->p,
                        ((struct delaunay_args *)args)->N,
                        ((struct delaunay_args *)args)->refs,
                        ((struct delaunay_args *)args)->dim,
                        ((struct delaunay_args *)args)->n_threads);
}
#else
DWORD delaunay_entry(void *args)
{
    ((struct delaunay_args *)args)->result = build_2dtree(
                        ((struct delaunay_args *)args)->p,
                        ((struct delaunay_args *)args)->N,
                        ((struct delaunay_args *)args)->refs,
                        ((struct delaunay_args *)args)->dim,
                        ((struct delaunay_args *)args)->n_threads);
    return 0;
}
#endif //__WIN32

gnu_attribute(malloc, nonnull)
struct delaunay_triangle_pool *triangulate(size_t N,
                                           const double nodes[static 2 * N],
                                           unsigned int n_threads)
{
    struct delaunay_triangle_pool *p = new_triangle_pool(2 * (N - 1));
    if (!p) {
        return NULL;
    }
    const double **refs = malloc(N * (sizeof *refs));
    if (!refs) {
        free_triangle_pool(p);
        return NULL;
    }
    for (size_t i = 0; i < N; i++) {
        refs[i] = nodes + 2 * i;
    }

    #ifdef __WIN32
    InitializeCriticalSection(&crit);
    #endif //__WIN32

    build_2dtree(p, N, refs, DIM_X, n_threads);
    free(refs);

    #ifdef __WIN32
    DeleteCriticalSection(&crit);
    #endif //__WIN32

    return p;
}
