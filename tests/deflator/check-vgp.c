#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mdwf.h>
#include "vfermion-test.h"

#define ESIZE 2
#define LS    3
#define UDIM  4

static int pattern[ESIZE * LS * UDIM * 10];

int
main(int argc, char *argv[])
{
    struct Fermion *f1 = new_fermion(ESIZE, LS);
    struct Fermion *f2 = new_fermion(ESIZE, LS);
    struct Fermion *f3 = new_fermion(ESIZE, LS);
    struct Fermion *ft = new_fermion(ESIZE, LS);
    struct vFermion *v0 = new_vfermion(ESIZE, LS, UDIM);
    struct vFermion *v1 = new_vfermion(ESIZE, LS, 2 * UDIM);
    int i;

    for (i = 0; i < sizeof (pattern) / sizeof (pattern[0]); i++) {
        pattern[i] = lrand48() % 18 - 9;
    }

    mk_fermion(f1, ESIZE,  LS, pattern, 24);
    show_fermion("F1", f1, ESIZE, LS);
    mk_fermion(f2, ESIZE,  LS, pattern, 20);
    show_fermion("F2", f2, ESIZE, LS);
    mk_fermion(f3, ESIZE,  LS, pattern, 13);
    show_fermion("F3", f3, ESIZE, LS);

    qx(vf_put)(ESIZE, LS, v0, UDIM, 0, f1);
    qx(vf_put)(ESIZE, LS, v0, UDIM, 1, f2);
    qx(vf_put)(ESIZE, LS, v0, UDIM, 2, f3);

    qx(vf_get)(ESIZE, LS, ft, v0, UDIM, 0);
    show_fermion("v0", ft, ESIZE, LS);

    qx(vf_get)(ESIZE, LS, ft, v0, UDIM, 2);
    show_fermion("v2", ft, ESIZE, LS);

    qx(vf_copy)(ESIZE, LS, 2, v1, 2 * UDIM, 0, v0, UDIM, 0);
    qx(vf_copy)(ESIZE, LS, 2, v1, 2 * UDIM, 2, v0, UDIM, 1);
    qx(vf_copy)(ESIZE, LS, 2, v1, 2 * UDIM, 4, v0, UDIM, 1);

    for (i = 0; i < UDIM; i++) {
        char n[10];

        qx(vf_get)(ESIZE, LS, ft, v0, UDIM, i);
        sprintf(n, "v0 %d", i);
        show_fermion(n, ft, ESIZE, LS);
    }

    for (i = 0; i < 2 * UDIM; i++) {
        char n[10];

        qx(vf_get)(ESIZE, LS, ft, v1, 2 * UDIM, i);
        sprintf(n, "zz %d", i);
        show_fermion(n, ft, ESIZE, LS);
    }

    return 0;
}
