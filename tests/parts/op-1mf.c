#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

int
op_1mF_odd(struct Fermion *r_odd,
	   struct Q(State) *state,
	   const struct SUn *U,
	   const struct Fermion *s_odd,
	   const struct Fermion *s_even)
{
    struct eo_lattice *xy = &state->odd;
    int Ls = state->Ls;

    if (q(setup_comm)(state, sizeof (REAL)))
	return 1;

    op_boundary(xy, Ls, up_project_n, down_project_n, U, s_even);

    if (xy->h_valid)
	QMP_start(xy->handle);

    qx(do_1mF)(r_odd, 0, xy->body_size, Ls,
	       xy->neighbor, U, s_odd, s_even, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    qx(do_1mF)(r_odd, xy->body_size, xy->face_size, Ls,
	       xy->neighbor, U, s_odd, s_even, xy->receive_buf);

    return 0;
}


int
op_1mF_even(struct Fermion *r_even,
	   struct Q(State) *state,
	   const struct SUn *U,
	   const struct Fermion *s_even,
	   const struct Fermion *s_odd)
{
    struct eo_lattice *xy = &state->even;
    int Ls = state->Ls;

    if (q(setup_comm)(state, sizeof (REAL)))
	return 1;

    op_boundary(xy, Ls, up_project_n, down_project_n, U, s_odd);

    if (xy->h_valid)
	QMP_start(xy->handle);

    qx(do_1mF)(r_even, 0, xy->body_size, Ls,
	       xy->neighbor, U, s_even, s_odd, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    qx(do_1mF)(r_even, xy->body_size, xy->face_size, Ls,
	       xy->neighbor, U, s_even, s_odd, xy->receive_buf);

    return 0;
}
