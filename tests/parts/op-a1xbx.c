#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

int
op_A1xBx_even(struct Fermion *r_e,
	      const struct Q(Parameters) *params,
	      const struct Fermion *s_e)
{
    qx(do_A1xBx)(r_e, params->state->even.full_size, params->state->even.Ls,
		 params->BxpTable,
		 params->BxmTable,
		 params->AxipTable,
		 params->AximTable,
		 s_e);
    return 0;
}

int
op_A1xBx_odd(struct Fermion *r_o,
	      const struct Q(Parameters) *params,
	      const struct Fermion *s_o)
{
    qx(do_A1xBx)(r_o, params->state->odd.full_size, params->state->odd.Ls,
		 params->BxpTable,
		 params->BxmTable,
		 params->AxipTable,
		 params->AximTable,
		 s_o);
    return 0;
}

int
op_A1xBx(struct QX(Fermion) *r,
	 const struct Q(Parameters) *params,
	 const struct QX(Fermion) *s)
{
    op_A1xBx_even(r->even, params, s->even);
    op_A1xBx_odd(r->odd, params, s->odd);
    return 0;
}
