#include <stdio.h>
#include <stdlib.h>
#define __USE_POSIX199309
#include <time.h>
#include "include/geometry.h"
#include "include/delaunay.h"

#define N 1234
#define RANDOMIZE_MIN -100.0
#define RANDOMIZE_MAX 100.0


#ifdef __WIN32
#define LONGF(s) "%ll"#s
#else
#define LONGF(s) "%l"#s
#endif //__WIN32


extern void write_triangulation(const struct delaunay_triangle_pool *p,
                                double xlow, double xhigh,
                                double ylow, double yhigh,
                                int width_px, int height_px,
                                const char *filename);

extern void write_voronoi(const struct delaunay_triangle_pool *p,
                          double xlow, double xhigh,
                          double ylow, double yhigh,
                          int width_px, int height_px,
                          const char *filename);


void randomize_node(double node[], double min, double max);

void print_time(const struct timespec *t0, const struct timespec *t1);


int main(int argc, char *argv[])
{
    unsigned int seed;
    if (argc > 1) {
        sscanf(argv[1], "%u", &seed);
    } else {
        seed = time(NULL);
    }
    printf("Seeding with %u\n", seed);
    srand(seed);

    double *nodes = malloc(2 * N * (sizeof *nodes));
    if (!nodes) {
        perror("Failed to allocate sites");
        return -1;
    }
    for (size_t i = 0; i < N; i++) {
        randomize_node(nodes + 2 * i, RANDOMIZE_MIN, RANDOMIZE_MAX);
    }
    /* double nodes[] = {
        -20.0, -20.0,
        20.0, -20.0,
        0.0, -20.0
    }; */

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    struct delaunay_triangle_pool *p = triangulate(N, nodes, 8);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (!p) {
        perror("Failed to allocate memory for triangulation");
        return -1;
    } else {
        puts("Triangulation complete!");
    }

    print_time(&t0, &t1);

    write_triangulation(p, -102, 102, -102, 102, 2048, 2048, "triangulation.jpg");
    write_voronoi(p, -102, 102, -102, 102, 2048, 2048, "voronoi.jpg");

    free_triangle_pool(p);
    free(nodes);

    #ifdef __WIN32
    system("pause");
    #endif //__WIN32
    
    return 0;
}


void randomize_node(double node[], double min, double max)
{
    for (size_t i = 0; i < 2; i++) {
        node[i] = ((double)rand() / RAND_MAX) * (max - min) + min;
    }
}

void print_time(const struct timespec *t0, const struct timespec *t1)
{
    time_t nsec = t1->tv_nsec - t0->tv_nsec;
    time_t sec = t1->tv_sec - t0->tv_sec;
    if (sec) {
        printf("Took %.3f s\n", ((double)sec * 1e9 + (double)nsec) / 1e9);
    } else {
        if (nsec < 1000) {
            printf("Took " LONGF(i) " ns\n", nsec);
        } else if (nsec < 1000000) {
            printf("Took " LONGF(i) " µs\n", nsec / 1000);
        } else if (nsec < 1000000000) {
            printf("Took %.3f ms\n", (double)nsec / 1e6);
        }
    }
}
