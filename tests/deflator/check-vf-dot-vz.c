#include <stdio.h>
#include <stdlib.h>
#include <mdwf.h>
#include <math.h>
#include "vfermion-test.h"

void
test_dot_zv(const char *name,
            int esize, 
            int ls,
            struct vFermion *v,
            int v_w, int v_b, int v_l,
            double *zv,
            struct Fermion *t0,
            struct Fermion *t1)
{
    int i;
    double f0;
    double ff;

    qx(vf_dot_vz)(esize, ls, t0, v, v_w, v_b, v_l, zv);
    qx(f_norm)(&f0, esize, ls, t0);
    for (i = 0; i < v_l; i++) {
        qx(vf_get)(esize, ls, t1, v, v_w, i + v_b);
        qx(f_cadd2)(t0, esize, ls, -zv[2*i], -zv[2*i+1], t1);
    }
    qx(f_norm)(&ff, esize, ls, t0);
    printf("DELTA: %15.7e\n", sqrt(ff/f0));
}

int
main(int argc, char *argv[])
{
    int esize;
    int ls;
    int width;
    int v_begin, v_len;
    struct vFermion *v0;
    struct Fermion *t0;
    struct Fermion *t1;

    if (argc != 6) {
        fprintf(stderr, "Usage: check/fv:zv esize ls width v_b v_l\n");
        return 1;
    }
    esize = atoi(argv[1]);
    ls = atoi(argv[2]);
    width = atoi(argv[3]);
    v_begin = atoi(argv[4]);
    v_len = atoi(argv[5]);

    t0 = new_fermion(esize, ls);
    t1 = new_fermion(esize, ls);
    v0 = new_vfermion(esize, ls, width);
    double zv[2 * width];
    construct_d(2 * width, zv, -756.345234);

    construct_vf(esize, ls, width, v0, t0, 12.45);
    construct_f(esize, ls, t0, 0.45673);

    test_dot_zv("ZV", esize, ls, v0, width, v_begin, v_len, zv, t0, t1);

    return 0;
}
