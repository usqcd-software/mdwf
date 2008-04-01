#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"
#include "op-routines.h"

char *op_a_name = "lib:DDW";
char *op_b_name = "conj(lib:DDW)";

double
read_gauge(int dir,
	   const int pos[4],
	   int a, int b,
	   int re_im,
	   void *env)
{
    int i;
    unsigned int v = sum_init(seed_u);

    v = sum_add(v, re_im);
    v = sum_add(v, a);
    v = sum_add(v, b);
    for (i = 0; i < 4; i++) {
	v = sum_add(v, pos[i]);
	v = sum_add(v, dir);
    }
    v = sum_add(v, seed_u);
    return sum_fini(v);
}

static double
read_fermion(unsigned int seed,
	     const int pos[5],
	     int c, int d,
	     int re_im)
{
    int i;
    unsigned int v = sum_init(seed);

    v = sum_add(v, c);
    v = sum_add(v, d);
    v = sum_add(v, re_im);
    for (i = 0; i < 5; i++)
	v = sum_add(v, pos[i]);
    v = sum_add(v, seed);
    return sum_fini(v);
}

double
read_fermion_a(const int pos[5],
	       int c, int d,
	       int re_im,
	       void *env)
{
    return read_fermion(seed_a, pos, c, d, re_im);
}

double
read_fermion_b(const int pos[5],
	       int c, int d,
	       int re_im,
	       void *env)
{
    return read_fermion(seed_b, pos, c, d, re_im);
}

int
operator_a(void)
{
    struct QX(Fermion) *fermion_x;
    double x, y;

    if (QOP_MDWF_allocate_fermion(&fermion_x, state)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QX(DDW_operator)(fermion_x, params, gauge, fermion_a)) {
	zprint("operator_DDWF(): failed");
	return 1;
    };

    dot_fermion(&x, &y, fermion_b, fermion_x);

    QOP_MDWF_free_fermion(&fermion_x);
    zprint("A+FB   : %20.10e %20.10e", x, y);
    return 0;
}

int
operator_b(void)
{
    struct QX(Fermion) *fermion_x;
    double x, y;

    if (QOP_MDWF_allocate_fermion(&fermion_x, state)) {
	zprint("operator_b(): alloc failed on x");
	return 1;
    }

    if (QX(DDW_operator_conjugated)(fermion_x, params, gauge, fermion_b)) {
	zprint("operator_DDWF_conjugate(): failed");
	return 1;
    };

    dot_fermion(&x, &y, fermion_x, fermion_a);

    QOP_MDWF_free_fermion(&fermion_x);
    zprint("A*+B*F*: %20.10e %20.10e", x, y);
    return 0;
}
