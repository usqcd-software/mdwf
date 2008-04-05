#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

int
op_BA1_odd(struct Fermion *r_odd,
	   const struct Q(Parameters) *params,
	   const struct Fermion *a_odd)
{
    struct Q(State) *state = params->state;
    struct eo_lattice *xy = &state->odd;
    int Ls = state->Ls;

    qx(do_BA1)(r_odd, xy->full_size, Ls,
		params->BpTable, params->BmTable,
		params->AipTable, params->AimTable,
		a_odd);
    return 0;
}

int
op_BA1_even(struct Fermion *r_even,
	    const struct Q(Parameters) *params,
	    const struct Fermion *a_even)
{
    struct Q(State) *state = params->state;
    struct eo_lattice *xy = &state->even;
    int Ls = state->Ls;

    qx(do_BA1)(r_even, xy->full_size, Ls,
		params->BpTable, params->BmTable,
		params->AipTable, params->AimTable,
		a_even);
    return 0;
}
