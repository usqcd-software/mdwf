#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

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

int
op_F_even(struct Fermion *result_even,
	  struct Q(State) *state,
	  const struct SUn *U,
	  const struct Fermion *src_odd)
{
    struct eo_lattice *even = &state->even;
    int Ls = state->Ls;
    int i;

    if (q(setup_comm)(state, sizeof (REAL)))
	return 1;

    for (i = 0; i < Q(DIM); i++) {
	if (even->send_up_size[i])
	    (up_project[i])(even->send_up_buf[i],
			    even->send_up_size[i], Ls,
			    even->up_pack[i], U, src_odd);
	if (even->send_down_size[i])
	    (down_project[i])(even->send_down_buf[i],
			      even->send_down_size[i], Ls,
			      even->down_pack[i], src_odd);
    }

    if (even->h_valid)
	QMP_start(even->handle);

    if (even->body_size)
	qx(do_F)(result_even, 0, even->body_size, Ls,
		 even->body_neighbor,
		 U, src_odd, NULL);

    if (even->h_valid)
	QMP_wait(even->handle);
    
    if (even->face_size)
	qx(do_F)(result_even, even->body_size, even->face_size, Ls,
		 even->body_neighbor,
		 U, src_odd, even->receive_buf);
    return 0;
}

int
op_F_odd(struct Fermion *result_odd,
	 struct Q(State) *state,
	 const struct SUn *U,
	 const struct Fermion *src_even)
{
    struct eo_lattice *odd = &state->odd;
    int Ls = state->Ls;
    int i;

    if (q(setup_comm)(state, sizeof (REAL)))
	return 1;

    for (i = 0; i < Q(DIM); i++) {
	if (odd->send_up_size[i])
	    (up_project[i])(odd->send_up_buf[i],
			    odd->send_up_size[i], Ls,
			    odd->up_pack[i], U, src_even);
	if (odd->send_down_size[i])
	    (down_project[i])(odd->send_down_buf[i],
			      odd->send_down_size[i], Ls,
			      odd->down_pack[i], src_even);
    }

    if (odd->h_valid)
	QMP_start(odd->handle);

    if (odd->body_size)
	qx(do_F)(result_odd, 0, odd->body_size, Ls,
		 odd->body_neighbor,
		 U, src_even, NULL);

    if (odd->h_valid)
	QMP_wait(odd->handle);
    
    if (odd->face_size)
	qx(do_F)(result_odd, odd->body_size, odd->face_size, Ls,
		 odd->body_neighbor,
		 U, src_even, odd->receive_buf);
    return 0;
}

int
op_F(struct Q(Fermion) *result,
     struct Q(State) *state,
     const struct Q(Gauge) *gauge,
     const struct Q(Fermion) *source)
{
    int status = 0;

    status |= op_F_even(result->even, state, gauge->data, source->odd);
    status |= op_F_odd(result->odd, state, gauge->data, source->odd);

    return status;
}
