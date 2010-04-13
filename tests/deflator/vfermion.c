#include <stdlib.h>
#include <math.h>
#include <mdwf.h>
#include "vfermion-test.h"

struct Fermion *
new_fermion(int size, int ls)
{
    void *ref = malloc(2 * qx(sizeof_fermion(size, ls)) + 127);
    return (void *)(((uintptr_t)ref + 127) & (~(uintptr_t)127));
}

struct vFermion *
new_vfermion(int size, int ls, int dim)
{
    void *ref = malloc(2 * qx(sizeof_vfermion(size, ls, dim)) + 127);
    return (void *)(((uintptr_t)ref + 127) & (~(uintptr_t)127));
}

void
mk_fermion(struct Fermion *f, int size, int ls, const int *data, int stride)
{
    int i, l, c, d;
    const int *p;
    double *q;
    double z[2 * QOP_MDWF_FERMION_DIM * QOP_MDWF_COLORS * ls];

    for (i = 0; i < size; i++, data += stride) {
        for (p = data, q = z, l = 0; l < ls; l++) {
            for (c = 0; c < QOP_MDWF_COLORS; c++) {
                for (d = 0; d < QOP_MDWF_FERMION_DIM; d++, p += 2, q += 2) {
                    q[0] = p[0];
                    q[1] = p[1];
                }
            }
        }
        qx(put_fermion)(f, i, ls, z);
    }
}

void
show_fermion_pt(const char *name, int k, int ls, const struct Fermion *f)
{
    double z[2 * QOP_MDWF_FERMION_DIM * QOP_MDWF_COLORS * ls];
    int l, c, d;

    qx(get_fermion)(z, f, k, ls);

    for (l = 0; l < ls; l++) {
        for (d = 0; d < QOP_MDWF_FERMION_DIM; d++) {
            printf("%s %d %d [%d]:", name, k, l, d);
            for (c = 0; c < QOP_MDWF_COLORS; c++) {
                int idx = c + QOP_MDWF_COLORS * (d + QOP_MDWF_FERMION_DIM * l);
                printf("  %8.3f %8.3f",
                       z[0 + 2 * idx],
                       z[1 + 2 * idx]);
            }
            printf("\n");
        }
    }
}

void
show_fermion(const char *name, const struct Fermion *f, int size, int ls)
{
    int i;

    for (i = 0; i < size; i++)
        show_fermion_pt(name, i, ls, f);
    printf("\n");
}

void
construct_f(int esize, int ls, struct Fermion *f, double m)
{
    int j, l, k;
    double z[2 * QOP_MDWF_COLORS * QOP_MDWF_FERMION_DIM * ls];

    for (j = 0; j < esize; j++) {
        for (l = 0; l < ls; l++) {
            for (k = 0; k < sizeof (z) / sizeof (z[0]); k++) {
#if 0
                z[k] = sin((k * m - (j + l) * (k + 2)) / (j + l + m));
#else
                z[k] = drand48();
#endif
            }
        }
        qx(put_fermion)(f, j, ls, z);
    }
}

void
construct_vf(int esize, int ls, int width,
             struct vFermion *dst,
             struct Fermion *t, double m)
{
    int i;

    for (i = 0; i < width; i++) {
        construct_f(esize, ls, t, cos((i * m)/(i + m + 4)));
        qx(vf_put)(esize, ls, dst, width, i, t);
    }
}

void
construct_d(int len,
            double *d,
            double m)
{
    int i;

    for (i = 0; i < len; i++) {
#if 0
        d[i] = sin((i + len)/(i * i + m + len));
#else
        d[i] = drand48();
#endif
    }
}
