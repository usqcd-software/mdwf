#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"

char *op_a_name = "1";
char *op_b_name = "inverse(A) * A";

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

static int
operator_A(struct QX(Fermion) *result,
	   const struct Q(Parameters) *params,
	   const struct QX(Gauge) *gauge,
	   const struct QX(Fermion) *fermion)
{
    qx(do_A)(result->even,
	     result->state->even.full_size,
	     result->state->even.Ls,
	     params->ApTable,
	     params->AmTable,
	     fermion->even);
    qx(do_A)(result->odd,
	     result->state->odd.full_size,
	     result->state->odd.Ls,
	     params->ApTable,
	     params->AmTable,
	     fermion->odd);
    return 0;
}

static int
operator_Ainv(struct QX(Fermion) *result,
	      const struct Q(Parameters) *params,
	      const struct QX(Gauge) *gauge,
	      const struct QX(Fermion) *fermion)
{
    qx(do_A_inverse)(result->even,
		     result->state->even.full_size,
		     result->state->even.Ls,
		     params->AipTable,
		     params->AimTable,
		     fermion->even);
    qx(do_A_inverse)(result->odd,
		     result->state->odd.full_size,
		     result->state->odd.Ls,
		     params->AipTable,
		     params->AimTable,
		     fermion->odd);
    return 0;
}

static void
dot(double *v_r, double *v_i,
    const struct QX(Fermion) *a,
    const struct QX(Fermion) *b)
{
    double r1, r2, i1, i2;

    qx(f_dot)(&r1, &i1, a->state->even.full_size, a->state->even.Ls,
	      a->even, b->even);
    qx(f_dot)(&r2, &i2, a->state->odd.full_size, a->state->odd.Ls,
	      a->odd, b->odd);
    *v_r = r1 + r2;
    *v_i = i1 + i2;
}

int
operator_a(void)
{
    double x, y;

    dot(&x, &y, fermion_b, fermion_a);

    zprint("dot     : %20.10e %20.10e", x, y);
    return 0;
}

int
operator_b(void)
{
    double x, y;
    struct QX(Fermion) *fermion_x;
    struct QX(Fermion) *fermion_y;

    if (QOP_MDWF_allocate_fermion(&fermion_x, gauge->state)) {
	zprint("operator_1/A A (): alloc 1 failed");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&fermion_y, gauge->state)) {
	zprint("operator_1/A A (): alloc 2 failed");
	return 1;
    }

    operator_A(fermion_x, params, gauge, fermion_a);
    operator_Ainv(fermion_y, params, gauge, fermion_x);
    dot(&x, &y, fermion_b, fermion_y);
    QOP_MDWF_free_fermion(&fermion_x);
    QOP_MDWF_free_fermion(&fermion_y);
    zprint("identity: %20.10e %20.10e", x, y);
    return 0;
}
