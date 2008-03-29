#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

int
op_A1_even(struct Fermion *result,
	   const struct Q(Parameters) *params,
	   const struct Fermion *fermion)
{
    qx(do_A_inverse)(result,
		     params->state->even.full_size,
		     params->state->even.Ls,
		     params->AipTable,
		     params->AimTable,
		     fermion);
    return 0;
}

int
op_A1_odd(struct Fermion *result,
	  const struct Q(Parameters) *params,
	  const struct Fermion *fermion)
{
    qx(do_A_inverse)(result,
		     params->state->odd.full_size,
		     params->state->odd.Ls,
		     params->AipTable,
		     params->AimTable,
		     fermion);
    return 0;
}

int
op_A1(struct QX(Fermion) *result,
      const struct Q(Parameters) *params,
      const struct QX(Fermion) *fermion)
{
    op_A1_even(result->even, params, fermion->even);
    op_A1_odd(result->odd, params, fermion->odd);
    return 0;
}
