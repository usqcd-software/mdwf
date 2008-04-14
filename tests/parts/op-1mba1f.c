#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

int
op_1mBA1F_odd(struct Fermion *r_odd,
	      struct Q(State) *state,
	      const struct Q(Parameters) *params,
	      const struct SUn *U,
	      const struct Fermion *a_odd,
	      const struct Fermion *a_even)
{
    struct eo_lattice *xy = &state->odd;
    int Ls = state->Ls;
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    (up_project_n[i])(xy->send_up_buf[i],
			      xy->send_up_size[i], Ls,
			      xy->up_pack[i], U, a_even);
	if (xy->send_down_size[i])
	    (down_project_n[i])(xy->send_down_buf[i],
				xy->send_down_size[i], Ls,
				xy->down_pack[i], a_even);
    }
    if (xy->h_valid)
	QMP_start(xy->handle);

    qx(do_1mBA1F)(r_odd, 0, xy->body_size, Ls,
		  params->BpTable, params->BmTable,
		  params->AipTable, params->AimTable,
		  xy->neighbor, U, a_odd, a_even, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    qx(do_1mBA1F)(r_odd, xy->body_size, xy->face_size, Ls,
		  params->BpTable, params->BmTable,
		  params->AipTable, params->AimTable,
		  xy->neighbor, U, a_odd, a_even, xy->receive_buf);
    return 0;
}
