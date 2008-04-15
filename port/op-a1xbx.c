#include <mdwf.h>

void
qx(op_A1xBx)(struct Fermion *r_x,
	     struct eo_lattice *xy,
	     const struct Q(Parameters) *params,
	     const struct Fermion *s_x,
	     long long *flops)
{
    *flops += qx(do_A1xBx)(r_x,
			   xy->full_size, xy->Ls,
			   params->BxpTable,
			   params->BxmTable,
			   params->AxipTable,
			   params->AximTable,
			   s_x);
}
