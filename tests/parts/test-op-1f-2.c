#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"
#include "op-routines.h"

char *op_a_name = "<o|lib:1mF*|e>";
char *op_b_name = "<o|1mF*|e>";

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
/***************************************************************************/
static Up_project up_project_x[Q(DIM)] = {
  qx(proj_Ucg0minus),
  qx(proj_Ucg1minus),
  qx(proj_Ucg2minus),
  qx(proj_Ucg3minus)
};

static Down_project down_project_x[Q(DIM)] = {
  qx(proj_g0plus),
  qx(proj_g1plus),
  qx(proj_g2plus),
  qx(proj_g3plus)
};

static void
compute_1mFx(struct Fermion *r_x,
	     struct Q(State) *state,
	     const struct SUn *U,
	     const struct Fermion *a_x,
	     const struct Fermion *a_y)
{
    struct eo_lattice *xy = &state->odd;
    int Ls = state->Ls;
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    (up_project_x[i])(xy->send_up_buf[i],
			      xy->send_up_size[i], Ls,
			      xy->up_pack[i], U, a_y);
	if (xy->send_down_size[i])
	    (down_project_x[i])(xy->send_down_buf[i],
				xy->send_down_size[i], Ls,
				xy->down_pack[i], a_y);
    }
    if (xy->h_valid)
	QMP_start(xy->handle);

    qx(do_1mFx)(r_x, 0, xy->body_size, Ls,
		xy->body_neighbor, U, a_x, a_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    qx(do_1mFx)(r_x, xy->body_size, xy->face_size, Ls,
		xy->face_neighbor, U, a_x, a_y, xy->receive_buf);
}

/***************************************************************************/

int
operator_a(void)
{
    struct QX(Fermion) *fermion_x;
    double x, y;

    if (QOP_MDWF_import_fermion(&fermion_x, state, read_fermion_a, NULL)) {
	zprint("operator_a(): alloc failed on x");
	return 1;
    }

    compute_1mFx(fermion_x->odd, state, gauge->data,
		 fermion_a->odd, fermion_a->even);

    dot_fermion(&x, &y, fermion_b, fermion_x);

    QOP_MDWF_free_fermion(&fermion_x);
    show("native", x, y);
    return 0;
}

/* same in parts */
static void
parts_1mfx(struct Fermion *r_o,
	   struct Q(State) *state,
	   const struct SUn *U,
	   struct Fermion *t0_o,
	   const struct Fermion *src_o,
	   const struct Fermion *src_e)
{
    op_Fx_odd(t0_o, state, U, src_e);
    qx(f_add3)(r_o, state->odd.full_size, state->odd.Ls, src_o, -1.0, t0_o);
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

    parts_1mfx(fermion_x->odd, state, gauge->data,
	       fermion_y->odd,
	       fermion_a->odd, fermion_a->even);

    dot_fermion(&x, &y, fermion_b, fermion_x);
    QOP_MDWF_free_fermion(&fermion_y);
    QOP_MDWF_free_fermion(&fermion_x);
    show("parts", x, y);
    return 0;
}
