#include <mdwf.h>

void
qx(op_A)(struct Fermion *r_x,
	 struct eo_lattice *xy,
	 const struct Q(Parameters) *params,
	 const struct Fermion *s_x,
	 long long *flops)
{
    *flops += qx(do_A)(r_x,
		       xy->full_size, xy->Ls,
		       params->ApTable,
		       params->AmTable,
		       s_x);
}
