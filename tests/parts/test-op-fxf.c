#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "opxtest.h"

char *op_a_name = "conj(F)";
char *op_b_name = "F";

/************************************************************************/
static Up_project up_project[Q(DIM)] = {
  qx(proj_Ucg0plus),
  qx(proj_Ucg1plus),
  qx(proj_Ucg2plus),
  qx(proj_Ucg3plus)
};

static Down_project down_project[Q(DIM)] = {
  qx(proj_g0minus),
  qx(proj_g1minus),
  qx(proj_g2minus),
  qx(proj_g3minus)
};

static void
compute_F(struct Q(State) *state,
	  struct eo_lattice *xy,
	  const struct Q(Parameters) *params,
	  struct Fermion *r_x,
	  const struct SUn *U,
	  const struct Fermion *s_y)
{
    int Ls = state->Ls;
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    (up_project[i])(xy->send_up_buf[i],
			    xy->send_up_size[i], Ls,
			    xy->up_pack[i], U, s_y);
	if (xy->send_down_size[i])
	    (down_project[i])(xy->send_down_buf[i],
			      xy->send_down_size[i], Ls,
			      xy->down_pack[i], s_y);
    }

    if (xy->h_valid)
	QMP_start(xy->handle);

    if (xy->body_size)
	qx(do_F)(r_x, 0, xy->body_size, Ls,
		 xy->body_neighbor,
		 U, s_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);
    
    if (xy->face_size)
	qx(do_F)(r_x, xy->body_size, xy->face_size, Ls,
		 xy->body_neighbor,
		 U, s_y, xy->receive_buf);
}

/************************************************************************/
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
compute_Fx(struct Q(State) *state,
	   struct eo_lattice *xy,
	   const struct Q(Parameters) *params,
	   struct Fermion *r_x,
	   const struct SUn *U,
	   const struct Fermion *s_y)
{
    int Ls = state->Ls;
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    (up_project_x[i])(xy->send_up_buf[i],
			      xy->send_up_size[i], Ls,
			      xy->up_pack[i], U, s_y);
	if (xy->send_down_size[i])
	    (down_project_x[i])(xy->send_down_buf[i],
				xy->send_down_size[i], Ls,
				xy->down_pack[i], s_y);
    }

    if (xy->h_valid)
	QMP_start(xy->handle);

    if (xy->body_size)
	qx(do_F_conj)(r_x, 0, xy->body_size, Ls,
		      xy->body_neighbor,
		      U, s_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);
    
    if (xy->face_size)
	qx(do_F_conj)(r_x, xy->body_size, xy->face_size, Ls,
		      xy->body_neighbor,
		      U, s_y, xy->receive_buf);
}

/*************************************************************************/
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
    struct Q(State) *state = result->state;

    if (q(setup_comm)(state, sizeof (REAL))) {
	zprint("setup_comm() failed");
	return 1;
    }
    compute_F(state, &state->even, params,
	      result->even, gauge->data, fermion->odd);
    compute_F(state, &state->odd, params,
	      result->odd, gauge->data, fermion->even);
    return 0;
}

static int
operator_Ax(struct QX(Fermion) *result,
	    const struct Q(Parameters) *params,
	    const struct QX(Gauge) *gauge,
	    const struct QX(Fermion) *fermion)
{
    struct Q(State) *state = result->state;

    if (q(setup_comm)(state, sizeof (REAL))) {
	zprint("setup_comm() failed");
	return 1;
    }
    compute_Fx(state, &state->even, params,
	       result->even, gauge->data, fermion->odd);
    compute_Fx(state, &state->odd, params,
	       result->odd, gauge->data, fermion->even);
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
    struct QX(Fermion) *fermion_x;

    if (QOP_MDWF_allocate_fermion(&fermion_x, gauge->state)) {
	zprint("operator_A(): alloc failed");
	return 1;
    }

    operator_A(fermion_x, params, gauge, fermion_a);
    dot(&x, &y, fermion_b, fermion_x);

    QOP_MDWF_free_fermion(&fermion_x);
    zprint("normal: %20.10e %20.10e", x, y);
    return 0;
}

int
operator_b(void)
{
    double x, y;
    struct QX(Fermion) *fermion_x;

    if (QOP_MDWF_allocate_fermion(&fermion_x, gauge->state)) {
	zprint("operator_Ax(): alloc failed");
	return 1;
    }

    operator_Ax(fermion_x, params, gauge, fermion_b);
    dot(&x, &y, fermion_x, fermion_a);
    QOP_MDWF_free_fermion(&fermion_x);
    zprint("conj  : %20.10e %20.10e", x, y);
    return 0;
}
