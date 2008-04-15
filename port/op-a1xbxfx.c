#include <mdwf.h>

void
qx(op_A1xBxFx)(struct Fermion *r_x,
	       struct eo_lattice *xy,
	       const struct Q(Parameters) *params,
	       const struct SUn *U,
	       const struct Fermion *s_y,
	       long long *flops,
	       long long *sent,
	       long long *received)
{
    int Ls = xy->Ls;

    qx(boundary)(xy, Ls, qx(up_project_x), qx(down_project_x), U, s_y, flops);
    
    if (xy->h_valid)
	QMP_start(xy->handle);

    *flops += qx(do_A1xBxFx)(r_x, 0, xy->body_size, Ls,
			     params->AxipTable, params->AximTable,
			     params->BxpTable, params->BxmTable,
			     xy->neighbor, U, s_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    *flops += qx(do_A1xBxFx)(r_x, xy->body_size, xy->face_size, Ls,
			     params->AxipTable, params->AximTable,
			     params->BxpTable, params->BxmTable,
			     xy->neighbor, U, s_y, xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;

}
