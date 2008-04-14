#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

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

    if (even->h_valid)
	QMP_start(even->handle);

    if (even->body_size)
	qx(do_F)(result_even, 0, even->body_size, Ls,
		 even->neighbor, U, src_odd, NULL);

    if (even->h_valid)
	QMP_wait(even->handle);

    if (even->face_size)
	qx(do_F)(result_even, even->body_size, even->face_size, Ls,
		 even->neighbor, U, src_odd, even->receive_buf);

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

    if (odd->h_valid)
	QMP_start(odd->handle);

    if (odd->body_size)
	qx(do_F)(result_odd, 0, odd->body_size, Ls,
		 odd->neighbor,
		 U, src_even, NULL);

    if (odd->h_valid)
	QMP_wait(odd->handle);

    if (odd->face_size)
	qx(do_F)(result_odd, odd->body_size, odd->face_size, Ls,
		 odd->neighbor,
		 U, src_even, odd->receive_buf);

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
