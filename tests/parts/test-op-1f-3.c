#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"
#include "op-routines.h"

char *op_a_name = "<|mlib:1mF|>";
char *op_b_name = "<|1mF*|>";

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

double
read_fermion_x(const int pos[5],
	       int c, int d,
	       int re_im,
	       void *env)
{
    return read_fermion(seed_a + seed_b, pos, c, d, re_im);
}

double
read_fermion_z(const int pos[5],
	       int c, int d,
	       int re_im,
	       void *env)
{
    return 0.0;
}

static void
show(const char *name, double x, double y)
{
    zprint("%-10s: %20.10e %20.10e", name, x, y);
}

static void
onemF(struct QX(Fermion) *r,
      const struct QX(Gauge) *U,
      const struct QX(Fermion) *t,
      const struct QX(Fermion) *s)
{
    struct Q(State) *state = U->state;
    long long flp = 0;
    long long snd = 0;
    long long rec = 0;

    qx(op_1mF)(r->even, &state->even, U->data, t->even, s->odd,
	       &flp, &snd, &rec);
    qx(op_1mF)(r->odd, &state->odd, U->data, t->odd, s->even,
	       &flp, &snd, &rec);
}

int
operator_a(void)
{
    struct QX(Fermion) *fermion_x;
    struct QX(Fermion) *fermion_y;
    double x, y;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_x, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QOP_MDWF_import_fermion(&fermion_y, state, read_fermion_z, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    onemF(fermion_y, gauge, fermion_x, fermion_a);

    dot_fermion(&x, &y, fermion_b, fermion_y);

    QOP_MDWF_free_fermion(&fermion_x);
    QOP_MDWF_free_fermion(&fermion_y);
    show("native", x, y);
    return 0;
}

/* same in parts */
static void
parts_1mf_o(struct Fermion *r_o,
	    struct Q(State) *state,
	    const struct SUn *U,
	    struct Fermion *t0_o,
	    const struct Fermion *src_o,
	    const struct Fermion *src_e)
{
    op_F_odd(t0_o, state, U, src_e);
    qx(f_add3)(r_o, state->odd.full_size, state->odd.Ls, src_o, -1.0, t0_o);
}

static void
parts_1mf_e(struct Fermion *r_e,
	    struct Q(State) *state,
	    const struct SUn *U,
	    struct Fermion *t0_e,
	    const struct Fermion *src_e,
	    const struct Fermion *src_o)
{
    op_F_even(t0_e, state, U, src_o);
    qx(f_add3)(r_e, state->even.full_size, state->even.Ls, src_e, -1.0, t0_e);
}

int
operator_b(void)
{
    double x, y;
    struct Q(State) *state = gauge->state;
    struct QX(Fermion) *fermion_x;
    struct QX(Fermion) *fermion_y;
    struct QX(Fermion) *fermion_z;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_x, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QOP_MDWF_import_fermion(&fermion_y, state, read_fermion_z, NULL)) {
	zprint("operator_a(): alloc failed on y");
	return 1;
    }

    if (QOP_MDWF_import_fermion(&fermion_z, state, read_fermion_z, NULL)) {
	zprint("operator_a(): alloc failed on y");
	return 1;
    }

    parts_1mf_o(fermion_y->odd, state, gauge->data,
		fermion_z->odd,
		fermion_x->odd, fermion_a->even);
    parts_1mf_e(fermion_y->even, state, gauge->data,
		fermion_z->even,
		fermion_x->even, fermion_a->odd);

    dot_fermion(&x, &y, fermion_b, fermion_y);
    QOP_MDWF_free_fermion(&fermion_z);
    QOP_MDWF_free_fermion(&fermion_y);
    QOP_MDWF_free_fermion(&fermion_x);
    show("parts", x, y);
    return 0;
}
