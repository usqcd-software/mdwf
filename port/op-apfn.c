#include <mdwf.h>

void
qx(op_ApF_norm)(struct Fermion *r_x,
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
    int Ls = xy->Ls;
    double n_b, n_f;

    qx(boundary)(xy, Ls, qx(up_project_n), qx(down_project_n), U, a_y, flops);

    if (xy->h_valid)
        QMP_start(xy->handle);

    *flops += qx(do_ApF_norm)(r_x, &n_b, 0, xy->body_size, Ls,
                              params->ApTable, params->AmTable,
                              xy->neighbor, U, a_x, a_y, NULL);

    if (xy->h_valid)
        QMP_wait(xy->handle);

    *flops += qx(do_ApF_norm)(r_x, &n_f, xy->body_size, xy->face_size, Ls,
                              params->ApTable, params->AmTable,
                              xy->neighbor, U, a_x, a_y, xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;
    *local_norm = n_b + n_f;
}

