#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"


static void
show_send(struct eo_lattice *xy)
{
    int i, k;
    double *v;
    
    xprint("boundary start");
    for (i = 0; i < 4; i++) {
	if (xy->send_up_size[i]) {
	    v = (double *)(xy->send_up_buf[i]);
	    xprint(" snd_up[%d] of %d", i, xy->send_up_size[i]);
	    for (k = 0; k < xy->Ls * 12 * xy->send_up_size[i]; k++, v++) {
		if (*v != 0.0)
		    xprint("  send_up[%4d][%4d]: %20.8f", i, k, *v);
	    }
	}
	if (xy->send_down_size[i]) {
	    v = (double *)(xy->send_down_buf[i]);
	    xprint(" snd_down[%d] of %d", i, xy->send_down_size[i]);
	    for (k = 0; k < xy->Ls * 12 * xy->send_down_size[i]; k++, v++) {
		if (*v != 0.0)
		    xprint("  send_down[%4d][%4d]: %20.8f", i, k, *v);
	    }
	}
    }
    xprint("boundary end");
}

static void
show_receive(struct eo_lattice *xy)
{
    int i, k;
    double *v;
    xprint("buffer start");
    for (i = 0; i < 4; i++) {
	if (xy->receive_up_size[i]) {
	    xprint(" rcv_up[%d] of %d", i, xy->receive_up_size[i]);
	    v = (double *)(xy->receive_buf[i]);
	    for (k = 0; k < xy->Ls * 12 * xy->receive_up_size[i]; k++, v++)
		if (*v != 0.0)
		    xprint("  rec_up[%d][%4d]: %20.8f", i, k, *v);
	}
	if (xy->receive_down_size[i]) {
	    xprint(" rcv_down[%d] of %d", i, xy->receive_down_size[i]);
	    v = (double *)(xy->receive_buf[i + 4]);
	    for (k = 0; k < xy->Ls * 12 * xy->receive_down_size[i]; k++, v++)
		if (*v != 0.0)
		    xprint("  rec_down[%d][%4d]: %20.8f", i, k, *v);
	}
    }
    xprint("buffer end");
}

int
op_F_even(struct Fermion *result_even,
	  struct Q(State) *state,
	  const struct SUn *U,
	  const struct Fermion *src_odd)
{
    struct eo_lattice *even = &state->even;
    int Ls = state->Ls;

    if (q(setup_comm)(state, sizeof (REAL)))
	return 1;

    op_boundary(even, Ls, up_project_n, down_project_n, U, src_odd);
    show_send(even);
    if (even->h_valid)
	QMP_start(even->handle);

    xprint("even: body %d face %d", even->body_size, even->face_size);
    if (even->body_size)
	qx(do_F)(result_even, 0, even->body_size, Ls,
		 even->body_neighbor,
		 U, src_odd, NULL);

    if (even->h_valid)
	QMP_wait(even->handle);

    show_receive(even);
/*
    receive_buf[0..3] - from up
    receive_buf[4..7] - from down
*/
#if 1
    if (even->face_size)
	qx(do_F)(result_even, even->body_size, even->face_size, Ls,
		 even->body_neighbor,
		 U, src_odd, even->receive_buf);
#endif
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

    if (q(setup_comm)(state, sizeof (REAL)))
	return 1;

    op_boundary(odd, Ls, up_project_n, down_project_n, U, src_even);
    show_send(odd);

    if (odd->h_valid)
	QMP_start(odd->handle);

    xprint("odd: body %d face %d", odd->body_size, odd->face_size);
    if (odd->body_size)
	qx(do_F)(result_odd, 0, odd->body_size, Ls,
		 odd->body_neighbor,
		 U, src_even, NULL);

    if (odd->h_valid)
	QMP_wait(odd->handle);

    show_receive(odd);

#if 1
    if (odd->face_size)
	qx(do_F)(result_odd, odd->body_size, odd->face_size, Ls,
		 odd->body_neighbor,
		 U, src_even, odd->receive_buf);
#endif
    return 0;
}

int
op_F(struct Q(Fermion) *result,
     const struct Q(Gauge) *gauge,
     const struct Q(Fermion) *source)
{
    struct Q(State) *state = gauge->state;
    int status = 0;

    status |= op_F_even(result->even, state, gauge->data, source->odd);
    status |= op_F_odd(result->odd, state, gauge->data, source->even);

    return status;
}
