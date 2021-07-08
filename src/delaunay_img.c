/**
 *  Quickly dumps the triangle pool to either a triangulation jpeg or a 
 *  Voronoi diagram. Uses the STB image library to write image files. This 
 *  file's only use is to make sure that the output is correct
**/
#include <stdio.h>
#include <stdlib.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "include/stb_image_write.h"
#include "include/geometry.h"
#include "include/delaunay.h"

#ifdef __GNUC__
#define gnu_attribute(...) __attribute__((__VA_ARGS__))
#else
#define gnu_attribute(...)
#endif


struct image_info {
    double min[2], max[2];
    double range[2];
    int res[2];
};

enum {
    DIM_X = 0,
    DIM_Y = 1
};

#define FEATHER 2 // expand the bounding box by this, in px

#define LINEWIDTH 0.02

gnu_attribute(nonnull)
/// Bounds has the following form:
/// ibegin, iend, istride
/// start iteration at ibegin to iend every level, increment by the x res 
/// each time, and do this istride times
///
static void find_bounding_box(const double *restrict v1, const double *restrict v2,
                              int bounds[static 4], const struct image_info *info)
{
    if (v1[DIM_X] < v2[DIM_X]) {
        bounds[0] = (v1[DIM_X] - info->min[DIM_X]) / (info->range[DIM_X]) * info->res[DIM_X];
        bounds[2] = (v2[DIM_X] - info->min[DIM_X]) / (info->range[DIM_X]) * info->res[DIM_X];
    } else {
        bounds[0] = (v2[DIM_X] - info->min[DIM_X]) / (info->range[DIM_X]) * info->res[DIM_X];
        bounds[2] = (v1[DIM_X] - info->min[DIM_X]) / (info->range[DIM_X]) * info->res[DIM_X];
    }
    bounds[0] -= FEATHER;
    if (bounds[0] < 0) {
        bounds[0] = 0;
    }
    bounds[2] += FEATHER;
    if (bounds[2] >= info->res[DIM_X]) {
        bounds[2] = info->res[DIM_X] - 1;
    }
    if (v1[DIM_Y] < v2[DIM_Y]) {
        bounds[1] = (v1[DIM_Y] - info->min[DIM_Y]) / info->range[DIM_Y] * info->res[DIM_Y];
        bounds[3] = (v2[DIM_Y] - info->min[DIM_Y]) / info->range[DIM_Y] * info->res[DIM_Y];
    } else {
        bounds[1] = (v2[DIM_Y] - info->min[DIM_Y]) / info->range[DIM_Y] * info->res[DIM_Y];
        bounds[3] = (v1[DIM_Y] - info->min[DIM_Y]) / info->range[DIM_Y] * info->res[DIM_Y];
    }
    bounds[1] -= FEATHER;
    if (bounds[1] < 0) {
        bounds[1] = 0;
    }
    bounds[3] += FEATHER;
    if (bounds[3] >= info->res[DIM_Y]) {
        bounds[3] = info->res[DIM_Y] - 1;
    }
}

gnu_attribute(nonnull)
static void get_pos_vector(double r[static restrict 2], const double src[static restrict 2],
                           int i, int j, const struct image_info *info)
{
    double scratch[2];
    scratch[DIM_X] = info->min[DIM_X] + (info->range[DIM_X] * (double)i / info->res[DIM_X]);
    scratch[DIM_Y] = info->min[DIM_Y] + (info->range[DIM_Y] * (double)j / info->res[DIM_Y]);
    relative_twovec(src, scratch, r);
}

gnu_attribute(nonnull)
static void draw_edge(float *raw, const double *restrict v1,
                      const double *restrict v2,
                      const struct image_info *info)
{
    int bounds[4]; // top left index, xrange, yrange
    find_bounding_box(v1, v2, bounds, info);
    double r[2];
    double b[2];
    relative_twovec(v1, v2, b);
    const double b_sqrd = dot_twovec(b, b);
    double proj, cross, d_sqrd;
    raw += bounds[1] * info->res[DIM_X];
    for (int j = bounds[1]; j <= bounds[3]; j++) {
        for (int i = bounds[0]; i <= bounds[2]; i++) {
            get_pos_vector(r, v1, i, j, info);
            proj = dot_twovec(r, b);
            if (proj <= 0) {
                d_sqrd = dot_twovec(r, r);
            } else if (proj >= b_sqrd) {
                r[DIM_X] -= v2[DIM_X];
                r[DIM_Y] -= v2[DIM_Y];
                d_sqrd = dot_twovec(r, r);
            } else {
                cross = cross_twovec(r, b);
                d_sqrd = (cross * cross) / b_sqrd;
            }
            if (d_sqrd < LINEWIDTH) {
                raw[i] += LINEWIDTH - d_sqrd;
            }
        }
        raw += info->res[DIM_X];
    }
}

gnu_attribute(nonnull, pure)
static float maximum_float(const float *raw, const struct image_info *info)
{
    float max = raw[0];
    for (int j = 0; j < info->res[DIM_Y]; j++) {
        for (int i = 0; i < info->res[DIM_X]; i++) {
            if (raw[i] > max) {
                max = raw[i];
            }
        }
        raw += info->res[DIM_X];
    }
    return max;
}

#define FADE 4

gnu_attribute(const)
static unsigned char norm_Y(unsigned char x)
{
    unsigned int y = x;
    for (int i = 0; i < FADE; i++) {
        y *= y;
        y /= 255;
    }
    if (y > 255) {
        y = 255;
    }
    return y;
}

gnu_attribute(nonnull)
static void normalize_image(const float *raw, unsigned char img[],
                            const struct image_info *info)
{
    float max = maximum_float(raw, info);
    raw += (info->res[DIM_Y] - 1) * info->res[DIM_X];
    for (int j = 0; j < info->res[DIM_Y]; j++) {
        for (int i = 0; i < info->res[DIM_X]; i++) {
            img[i] = norm_Y(255.0f - 255.0f * raw[i] / max);
        }
        raw -= info->res[DIM_X];
        img += info->res[DIM_X];
    }
}

gnu_attribute(nonnull, malloc)
static unsigned char *create_delaunay(const struct delaunay_triangle_pool *p,
                                      const struct image_info *info)
{
    size_t N_px = info->res[DIM_Y] * info->res[DIM_X];
    float *raw = calloc(N_px, sizeof *raw);
    if (!raw) {
        return NULL;
    }
    for (const struct delaunay_triangle *t = p->start; t < p->end; t++) {
        draw_edge(raw, t->vertices[0], t->vertices[1], info);
        if (t->vertices[2]) {
            draw_edge(raw, t->vertices[1], t->vertices[2], info);
            draw_edge(raw, t->vertices[2], t->vertices[0], info);
        }
    }
    unsigned char *img = calloc(N_px, sizeof *img);
    if (!img) {
        goto create_image_end;
    }
    normalize_image(raw, img, info);
create_image_end:
    free(raw);
    return img;
}

gnu_attribute(nonnull, malloc)
static unsigned char *create_voronoi(const struct delaunay_triangle_pool *p,
                                     const struct image_info *info)
{
    float *raw = calloc(info->res[DIM_Y] * info->res[DIM_X], sizeof *raw);
    if (!raw) {
        return NULL;
    }
    for (const struct delaunay_triangle *t = p->start; t < p->end; t++) {
        if (t->vertices[2]) {
            for (int i = 0; i < 3; i++) {
                if (t->adjacents[i].next->vertices[2]) {
                    draw_edge(raw, t->circumcenter, t->adjacents[i].next->circumcenter, info);
                }
            }
        }
    }
    unsigned char *img = calloc(info->res[DIM_X] * info->res[DIM_Y], sizeof *img);
    if (!img) {
        goto create_image_end;
    }
    normalize_image(raw, img, info);
create_image_end:
    free(raw);
    return img;
}

gnu_attribute(nonnull)
void write_triangulation(const struct delaunay_triangle_pool *p,
                         double xlow, double xhigh, double ylow, double yhigh,
                         int res_x, int res_y, const char *filename)
{
    struct image_info info;
    info.min[DIM_X] = xlow;
    info.max[DIM_X] = xhigh;
    info.range[DIM_X] = xhigh - xlow;
    info.res[DIM_X] = res_x;
    info.min[DIM_Y] = ylow;
    info.max[DIM_Y] = yhigh;
    info.range[DIM_Y] = yhigh - ylow;
    info.res[DIM_Y] = res_y;
    unsigned char *img = create_delaunay(p, &info);
    stbi_write_jpg(filename, res_x, res_y, 1, img, 90);
    free(img);
}

gnu_attribute(nonnull)
void write_voronoi(const struct delaunay_triangle_pool *p,
                   double xlow, double xhigh, double ylow, double yhigh,
                   int res_x, int res_y, const char *filename)
{
    struct image_info info;
    info.min[DIM_X] = xlow;
    info.max[DIM_X] = xhigh;
    info.range[DIM_X] = xhigh - xlow;
    info.res[DIM_X] = res_x;
    info.min[DIM_Y] = ylow;
    info.max[DIM_Y] = yhigh;
    info.range[DIM_Y] = yhigh - ylow;
    info.res[DIM_Y] = res_y;
    unsigned char *img = create_voronoi(p, &info);
    stbi_write_jpg(filename, res_x, res_y, 1, img, 90);
    free(img);
}
