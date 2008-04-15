#include <mdwf.h>
void
qx(op_1mFx)(struct Fermion *r_x,
	    struct eo_lattice *xy,
	    const struct SUn *U,
	    const struct Fermion *a_x,
	    const struct Fermion *a_y,
	    long long *flops,
	    long long *sent,
	    long long *received)
{
    int Ls = xy->Ls;

    qx(boundary)(xy, Ls, qx(up_project_x), qx(down_project_x), U, a_y, flops);

    if (xy->h_valid)
	QMP_start(xy->handle);

    *flops += qx(do_1mFx)(r_x, 0, xy->body_size, Ls,
			  xy->neighbor, U, a_x, a_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    *flops += qx(do_1mFx)(r_x, xy->body_size, xy->face_size, Ls,
			  xy->neighbor, U, a_x, a_y, xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;

}
