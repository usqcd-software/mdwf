#include <mdwf.h>

void
qx(op_BA1)(struct Fermion *r_x,
	   struct eo_lattice *xy,
	   const struct Q(Parameters) *params,
	   const struct Fermion *s_x,
	   long long *flops)
{
    int Ls = xy->Ls;

    *flops += qx(do_BA1)(r_x, xy->full_size, Ls,
			 params->BpTable, params->BmTable,
			 params->AipTable, params->AimTable,
			 s_x);
}
