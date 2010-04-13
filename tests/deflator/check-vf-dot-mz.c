#include <stdio.h>
#include <stdlib.h>
#include <mdwf.h>
#include <math.h>
#include "vfermion-test.h"

void
test_dot_zm(const char *name,
            int esize,
            int ls,
            struct vFermion *w,
            int w_w, int w_b, int w_l,
            struct vFermion *v,
            int v_w, int v_b, int v_l,
            double *z, int ldz,
            struct Fermion *t0,
            struct Fermion *t1)
{
    int i, j;
    double f0;
    double ff;
    double err;
    double m_e = 0;

    qx(vf_dot_mz)(esize, ls,
                  w, w_w, w_b, w_l,
                  v, v_w, v_b, v_l,
                  z, ldz);
    for (i = 0; i < w_l; i++) {
        qx(vf_get)(esize, ls, t0, w, w_w, i + w_b);
        qx(f_norm)(&f0, esize, ls, t0);
        for (j = 0; j < v_l; j++) {
            qx(vf_get)(esize, ls, t1, v, v_w, j + v_b);
            qx(f_cadd2)(t0, esize, ls,
                        -z[2*(j + ldz * i)], -z[2*(j+ ldz * i) + 1],
                        t1);
        }
        qx(f_norm)(&ff, esize, ls, t0);
        err = sqrt(ff/f0);
        printf("%s %3d %15.7e %15.7e %15.7e\n", name, i, err, ff, f0);
        if (m_e < err)
            m_e = err;
    }
    printf("MAX err: %15.7e\n", err);
}

int
main(int argc, char *argv[])
{
    int esize;
    int ls;
    int v_w, v_b, v_l;
    int w_w, w_b, w_l;
    struct vFermion *Vv;
    struct vFermion *Vw;
    struct Fermion *t0;
    struct Fermion *t1;
    int ldz;

    if (argc != 10) {
        fprintf(stderr, 
                "Usage: check/fv:zm esize ls v_w v_b v_l w_w w_b w_l ldz\n");
        return 1;
    }
    esize = atoi(argv[1]);
    ls = atoi(argv[2]);
    v_w = atoi(argv[3]);
    v_b = atoi(argv[4]);
    v_l = atoi(argv[5]);
    w_w = atoi(argv[6]);
    w_b = atoi(argv[7]);
    w_l = atoi(argv[8]);
    ldz = atoi(argv[9]);

    t0 = new_fermion(esize, ls);
    t1 = new_fermion(esize, ls);
    Vv = new_vfermion(esize, ls, v_w);
    Vw = new_vfermion(esize, ls, w_w);
    double z[2 * ldz * v_l];
    construct_d(2 * ldz * v_l, z, -756.345234);
    construct_vf(esize, ls, w_w, Vw, t0, 12.45);

    test_dot_zm("ZM", esize, ls,
                Vv, v_w, v_b, v_l,
                Vw, w_w, w_b, w_l,
                z, ldz, t0, t1);

    return 0;
}
