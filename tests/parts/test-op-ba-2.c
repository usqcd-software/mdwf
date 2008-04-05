#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"
#include "op-routines.h"

char *op_a_name = "<o|lib:B 1/A|e>";
char *op_b_name = "<o|B 1/A|e>";

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
    double x, y;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_a, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    op_BA1_odd(fermion_x->odd, params, fermion_a->odd);
    op_BA1_even(fermion_x->even, params, fermion_a->even);

    dot_fermion(&x, &y, fermion_b, fermion_x);

    QOP_MDWF_free_fermion(&fermion_x);
    show("native", x, y);
    return 0;
}

/* same in parts */
static void
parts_ba1_odd(struct Fermion *r_odd,
	      struct Q(Parameters) *params,
	      struct Fermion *tmp_odd,
	      const struct Fermion *src_odd)
{
    op_A1_odd(tmp_odd, params, src_odd);
    op_B_odd(r_odd, params, tmp_odd);
}

static void
parts_ba1_even(struct Fermion *r_even,
	       struct Q(Parameters) *params,
	       struct Fermion *tmp_even,
	       const struct Fermion *src_even)
{
    op_A1_even(tmp_even, params, src_even);
    op_B_even(r_even, params, tmp_even);
}

int
operator_b(void)
{
    double x, y;
    struct Q(State) *state = gauge->state;
    struct QX(Fermion) *fermion_x;
    struct QX(Fermion) *fermion_y;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_a, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QOP_MDWF_allocate_fermion(&fermion_y, state)) {
	zprint("operator_b(): alloc failed on y");
	return 1;
    }

    parts_ba1_odd(fermion_x->odd, params, fermion_y->odd, fermion_a->odd);
    parts_ba1_even(fermion_x->even, params, fermion_y->even, fermion_a->even);

    dot_fermion(&x, &y, fermion_b, fermion_x);
    QOP_MDWF_free_fermion(&fermion_y);
    QOP_MDWF_free_fermion(&fermion_x);
    show("parts", x, y);
    return 0;
}
