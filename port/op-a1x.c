#include <mdwf.h>

void
qx(op_A1x)(struct Fermion *r_x,
	   struct eo_lattice *xy,
	   const struct Q(Parameters) *params,
	   const struct Fermion *s_x,
	   long long *flops)
{
    *flops += qx(do_A_conj_inverse)(r_x,
				    xy->full_size, xy->Ls,
				    params->AxipTable,
				    params->AximTable,
				    s_x);
}
