#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"
#include "op-routines.h"

char *op_a_name = "<o|mlib:M|e>";
char *op_b_name = "<o|M|e>";

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

static double
read_fermion_x(const int pos[5],
	       int c, int d,
	       int re_im,
	       void *env)
{
    int *seed = env;
    return read_fermion(*seed, pos, c, d, re_im);
}

static void
show(const char *name, double x, double y)
{
    zprint("%-10s: %20.10e %20.10e", name, x, y);
}

int
operator_a(void)
{
    struct QX(HalfFermion) *f_a;
    struct QX(HalfFermion) *f_b;
    struct QX(HalfFermion) *f_x;
    double x, y;

    if (QX(import_half_fermion)(&f_a, state, read_fermion_x, &seed_a)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QX(import_half_fermion)(&f_b, state, read_fermion_x, &seed_b)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QX(import_half_fermion)(&f_x, state, read_fermion_x, &seed_a)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QX(M_operator)(f_x, params, gauge, f_a)) {
	zprint("operator_a(): M_operator failed");
	return 1;
    }

    QX(dot_half_fermion)(&x, &y, f_b, f_x);

    QOP_MDWF_free_half_fermion(&f_x);
    QOP_MDWF_free_half_fermion(&f_a);
    QOP_MDWF_free_half_fermion(&f_b);
    show("native", x, y);
    return 0;
}

/* same in parts */
static void
parts_m(struct Fermion *r_e,
	struct Q(State) *state,
	struct Q(Parameters) *params,
	const struct SUn *U,
	const struct Fermion *s_e,
	struct Fermion *t0_e,
	struct Fermion *t1_e,
	struct Fermion *t0_o,
	struct Fermion *t1_o)
	
{
    long long flops = 0;
    long long sent = 0;
    long long received = 0;

    qx(op_F)(t0_o, &state->odd, U, s_e, &flops, &sent, &received);
    qx(op_A1)(t1_o, &state->odd, params, t0_o, &flops);
    qx(op_B)(t0_o, &state->odd, params, t1_o, &flops);
    qx(op_F)(t0_e, &state->even, U, t0_o, &flops, &sent, &received);
    qx(op_A1)(t1_e, &state->even, params, t0_e, &flops);
    qx(op_B)(t0_e, &state->even, params, t1_e, &flops);
    qx(f_add3)(r_e, state->even.full_size, state->odd.Ls, s_e, -1.0, t0_e);
}

int
operator_b(void)
{
    double x, y;
    struct Q(State) *state = gauge->state;
    struct QX(HalfFermion) *f_a;
    struct QX(HalfFermion) *f_b;
    struct QX(HalfFermion) *f_x;
    struct QX(Fermion) *f_y;
    struct QX(Fermion) *f_z;


    if (QX(import_half_fermion)(&f_a, state, read_fermion_x, &seed_a)) {
	zprint("operator_b(): alloc failed on a");
	return 1;
    }

    if (QX(import_half_fermion)(&f_b, state, read_fermion_x, &seed_b)) {
	zprint("operator_b(): alloc failed on b");
	return 1;
    }

    if (QOP_MDWF_allocate_half_fermion(&f_x, state)) {
	zprint("operator_b(): alloc failed on x");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&f_y, state)) {
	zprint("operator_b(): alloc failed on y");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&f_z, state)) {
	zprint("operator_b(): alloc failed on z");
	return 1;
    }

    parts_m(f_x->even, state, params, gauge->data, f_a->even,
	    f_y->even, f_z->even, f_y->odd, f_z->odd);

    QX(dot_half_fermion)(&x, &y, f_b, f_x);
    QOP_MDWF_free_fermion(&f_z);
    QOP_MDWF_free_fermion(&f_y);
    QOP_MDWF_free_half_fermion(&f_x);
    QOP_MDWF_free_half_fermion(&f_a);
    QOP_MDWF_free_half_fermion(&f_b);
    show("parts", x, y);
    return 0;
}
