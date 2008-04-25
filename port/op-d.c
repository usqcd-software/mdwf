#include <mdwf.h>

void
qx(op_D)(struct Fermion *r_x,
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
    qx(op_ApF)(r_x, xy, params, U, a_x, tmp_y, flops, sent, received);
}
