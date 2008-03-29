#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

int
op_Bx_even(struct Fermion *result,
	   const struct Q(Parameters) *params,
	   const struct Fermion *fermion)
{
    qx(do_A)(result,
	     params->state->even.full_size,
	     params->state->even.Ls,
	     params->BxTable,
	     fermion);
    return 0;
}

int
op_Bx_odd(struct Fermion *result,
	  const struct Q(Parameters) *params,
	  const struct Fermion *fermion)
{
    qx(do_A)(result,
	     params->state->odd.full_size,
	     params->state->odd.Ls,
	     params->BxTable,
	     fermion);
    return 0;
}

int
op_Bx(struct QX(Fermion) *result,
      const struct Q(Parameters) *params,
      const struct QX(Fermion) *fermion)
{
    op_Bx_even(result->even, params, fermion->even);
    op_Bx_odd(result->odd, params, fermion->odd);
    return 0;
}
