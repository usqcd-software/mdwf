#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

int
op_1mFx_odd(struct Fermion *r_odd,
	    struct Q(State) *state,
	    const struct SUn *U,
	    const struct Fermion *a_odd,
	    const struct Fermion *a_even)
{
    struct eo_lattice *xy = &state->odd;
    int Ls = state->Ls;

    if (q(setup_comm)(state, sizeof (REAL)))
	return 1;

    op_boundary(xy, Ls, up_project_x, down_project_x, U, a_even);

    if (xy->h_valid)
	QMP_start(xy->handle);

    qx(do_1mFx)(r_odd, 0, xy->body_size, Ls,
		xy->neighbor, U, a_odd, a_even, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    qx(do_1mFx)(r_odd, xy->body_size, xy->face_size, Ls,
		xy->neighbor, U, a_odd, a_even, xy->receive_buf);

    return 0;
}
