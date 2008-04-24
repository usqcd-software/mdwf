#include <mdwf.h>

void
qx(op_1mBA1F_norm)(struct Fermion *r_x,
		   double *local_norm,
		   struct eo_lattice *xy,
		   const struct Q(Parameters) *params,
		   const struct SUn *U,
		   const struct Fermion *a_x,
		   const struct Fermion *a_y,
		   long long *flops,
		   long long *sent,
		   long long *received)
{
    double n_body, n_face;
    int Ls = xy->Ls;

    qx(boundary)(xy, Ls, qx(up_project_n), qx(down_project_n), U, a_y, flops);

    if (xy->h_valid)
	QMP_start(xy->handle);

    *flops += qx(do_1mBA1F_norm)(r_x, &n_body, 0, xy->body_size, Ls,
				 params->BpTable, params->BmTable,
				 params->AipTable, params->AimTable,
				 xy->neighbor, U, a_x, a_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    *flops += qx(do_1mBA1F_norm)(r_x, &n_face, xy->body_size, xy->face_size, Ls,
				 params->BpTable, params->BmTable,
				 params->AipTable, params->AimTable,
				 xy->neighbor, U, a_x, a_y, xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;
    *local_norm = n_body + n_face;
}
