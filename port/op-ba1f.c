#include <mdwf.h>

void
qx(op_BA1F)(struct Fermion *r_x,
	    struct eo_lattice *xy,
	    const struct Q(Parameters) *params,
	    const struct SUn *U,
	    const struct Fermion *s_y,
	    long long *flops,
	    long long *sent,
	    long long *received)
{
    int Ls = xy->Ls;

    qx(boundary)(xy, Ls, qx(up_project_n), qx(down_project_n), U, s_y, flops);

    if (xy->h_valid)
	QMP_start(xy->handle);

    *flops += qx(do_BA1F)(r_x, 0, xy->body_size, Ls,
			  params->BpTable, params->BmTable,
			  params->AipTable, params->AimTable,
			  xy->neighbor, U, s_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    *flops += qx(do_BA1F)(r_x, xy->body_size, xy->face_size, Ls,
			  params->BpTable, params->BmTable,
			  params->AipTable, params->AimTable,
			  xy->neighbor, U, s_y, xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;
}
