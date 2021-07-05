#pragma once

#ifndef DELAUNAY_H
#define DELAUNAY_H


/**
 *  A complete triangulation will only ever have 2 * n - 2 triangles. Since this 
 *  algorithm never actually deletes them, I use a stack allocator here.
 * 
 *  Consider a triangulation embedded on the surface of a closed manifold. If we
 *  include manifold triangles ("ghost" triangles) then each edge is bordered by
 *  two triangles, guaranteed. Since each triangle is guaranteed
 *  to contain three edges:
 *                                 2e = 3t.
 * 
 *  A triangulation on a closed manifold forms a spherical polyhedron. The
 *  Euler characteristic
 *                          (n + 1) - e + t = 2
 * 
 *  shows this triangle limit immediately (Dwyer, R. A. (1986)). In Dwyer's
 *  paper, he added 1 to t to account for the manifold face. In this case, I
 *  add 1 to n instead, to account for the manifold vertex.
**/


struct delaunay_triangle_pool {
    struct delaunay_triangle *first;
    struct delaunay_triangle {
        const double *vertices[3];
        double circumcenter[3];
        struct {
            struct delaunay_triangle *next;
            int next_side;
        } adjacents[3];
    } triangles[];
};


struct delaunay_triangle_pool *triangulate(size_t N,
                                           const double nodes[static 2 * N],
                                           unsigned int n_threads);

void free_triangle_pool(struct delaunay_triangle_pool *p);


#endif //DELAUNAY_H