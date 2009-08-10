#include <mdwf.h>

void
qx(op_D_norm)(struct Fermion *r_x,
              double *local_norm,
              struct eo_lattice *xy,
              struct eo_lattice *yx,
              const struct Q(Parameters) *params,
              const struct SUn *U,
              const struct Fermion *a_x,
              const struct Fermion *a_y,
              long long *flops,
              long long *sent,
              long long *received,
              struct Fermion *tmp_y)
{
    qx(op_B)(tmp_y, yx, params, a_y, flops);
    qx(op_ApF_norm)(r_x, local_norm,
                    xy, params, U, a_x, tmp_y, flops, sent, received);
}
