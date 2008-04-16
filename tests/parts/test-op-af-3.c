#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"
#include "op-routines.h"

char *op_a_name = "<o|mlib:ApF|e>";
char *op_b_name = "<o|ApF|e>";

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

static void
show(const char *name, double x, double y)
{
    zprint("%-10s: %20.10e %20.10e", name, x, y);
}

int
operator_a(void)
{
    struct QX(Fermion) *fermion_x;
    long long flt = 0;
    long long snd = 0;
    long long rec = 0;
    double x, y;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_a, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    qx(op_ApF)(fermion_x->odd, &state->odd, params, gauge->data,
	       fermion_a->odd, fermion_a->even, &flt, &snd, &rec);

    dot_fermion(&x, &y, fermion_b, fermion_x);

    QOP_MDWF_free_fermion(&fermion_x);
    show("native", x, y);
    return 0;
}

/* same in parts */
static void
parts_ba1f(struct Fermion *r_o,
	   struct Q(State) *state,
	   struct Q(Parameters) *params,
	   const struct SUn *U,
	   struct Fermion *t0_o, struct Fermion *t1_o,
	   const struct Fermion *src_o,
	   const struct Fermion *src_e)
{
    op_F_odd(t0_o, state, U, src_e);
    op_A_odd(t1_o, params, src_o);
    qx(f_add3)(r_o, state->odd.full_size, state->odd.Ls,
	       t0_o, +1.0, t1_o);
}

int
operator_b(void)
{
    double x, y;
    struct Q(State) *state = gauge->state;
    struct QX(Fermion) *fermion_x;
    struct QX(Fermion) *fermion_y;
    struct QX(Fermion) *fermion_z;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_a, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&fermion_y, state)) {
	zprint("operator_b(): alloc failed on y");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&fermion_z, state)) {
	zprint("operator_b(): alloc failed on z");
	return 1;
    }

    parts_ba1f(fermion_x->odd, state, params, gauge->data,
	       fermion_y->odd, fermion_z->odd,
	       fermion_a->odd, fermion_a->even);

    dot_fermion(&x, &y, fermion_b, fermion_x);
    QOP_MDWF_free_fermion(&fermion_z);
    QOP_MDWF_free_fermion(&fermion_y);
    QOP_MDWF_free_fermion(&fermion_x);
    show("parts", x, y);
    return 0;
}
