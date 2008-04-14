#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

int
op_BA1F_odd(struct Fermion *r_odd,
	    struct Q(State) *state,
	    const struct Q(Parameters) *params,
	    const struct SUn *U,
	    const struct Fermion *a_even)
{
    struct eo_lattice *xy = &state->odd;
    int Ls = state->Ls;

    if (q(setup_comm)(state, sizeof (REAL)))
	return 1;

    op_boundary(xy, Ls, up_project_n, down_project_n, U, a_even);

    if (xy->h_valid)
	QMP_start(xy->handle);

    qx(do_BA1F)(r_odd, 0, xy->body_size, Ls,
		params->BpTable, params->BmTable,
		params->AipTable, params->AimTable,
		xy->neighbor, U, a_even, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    qx(do_BA1F)(r_odd, xy->body_size, xy->face_size, Ls,
		params->BpTable, params->BmTable,
		params->AipTable, params->AimTable,
		xy->neighbor, U, a_even, xy->receive_buf);

    return 0;
}
