#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "optest.h"

char *op_name = "Operator F0";

/* Conservative version with contained communications */
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


static int
operator_F(struct QX(Fermion) *result,
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

double
read_gauge(int dir,
	   const int pos[4],
	   int a, int b,
	   int re_im,
	   void *env)
{
    return 0.0;
}

double
read_fermion(const int pos[5],
	     int c, int d,
	     int re_im,
	     void *env)
{
    int i;
    for (i = 0; i < 5; i++) {
	if (pos[i] != fermion_pos[i])
	    return 0.0;
    }
    if (c != fermion_color ||
	d != fermion_dirac ||
	re_im != fermion_reim)
	return 0.0;
    return 1.0;
}

void
write_fermion(const int pos[5],
	      int c, int d,
	      int re_im,
	      double value,
	      void *env)
{
    if (value != 0.0)
	xprint("   fermion[%d,%d,%d,%d;%d][%d,%d].%d = %g",
	       pos[0], pos[1], pos[2], pos[3], pos[4],
	       c, d, re_im, value);
}

int operator(void)
{
    return operator_F(result, params, gauge, fermion);
}
