#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"
#include "op-routines.h"

char *op_a_name = "<o|lib:A1*B*|e>";
char *op_b_name = "<o|A1*B*|e>";

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
/***************************************************************************/
static void
compute_A1xBx(struct QX(Fermion) *r,
	      struct Q(State) *state,
	      const struct Q(Parameters) *params,
	      const struct QX(Fermion) *s)
{
    qx(do_A1xBx)(r->even, state->even.full_size, state->even.Ls,
		 params->BxpTable,
		 params->BxmTable,
		 params->AxipTable,
		 params->AximTable,
		 s->even);
    qx(do_A1xBx)(r->odd, state->odd.full_size, state->odd.Ls,
		 params->BxpTable,
		 params->BxmTable,
		 params->AxipTable,
		 params->AximTable,
		 s->odd);
}
/***************************************************************************/

int
operator_a(void)
{
    struct QX(Fermion) *fermion_x;
    double x, y;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_z, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    compute_A1xBx(fermion_x, state, params, fermion_a);

    dot_fermion(&x, &y, fermion_b, fermion_x);

    QOP_MDWF_free_fermion(&fermion_x);
    show("native", x, y);
    return 0;
}

/* same in parts */
static void
parts_a1b(struct QX(Fermion) *r,
	  struct Q(State) *state,
	  struct Q(Parameters) *params,
	  struct QX(Fermion) *t,
	  const struct QX(Fermion) *src)
{
    op_Bx(t, params, src);
    op_A1x(r, params, t);
}

int
operator_b(void)
{
    double x, y;
    struct Q(State) *state = gauge->state;
    struct QX(Fermion) *fermion_x;
    struct QX(Fermion) *fermion_y;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_z, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    if (QOP_MDWF_import_fermion(&fermion_y, state, read_fermion_z, NULL)) {
	zprint("operator_a(): alloc failed on y");
	return 1;
    }

    parts_a1b(fermion_x, state, params, fermion_y, fermion_a);

    dot_fermion(&x, &y, fermion_b, fermion_x);
    QOP_MDWF_free_fermion(&fermion_x);
    QOP_MDWF_free_fermion(&fermion_y);
    show("parts", x, y);
    return 0;
}
