#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

void
dot_fermion(double *v_r, double *v_i,
	    const struct QX(Fermion) *a,
	    const struct QX(Fermion) *b)
{
    double r1, r2, i1, i2;

    qx(f_dot)(&r1, &i1, a->state->even.full_size, a->state->even.Ls,
		    a->even, b->even);
    qx(f_dot)(&r2, &i2, a->state->odd.full_size, a->state->odd.Ls,
		    a->odd, b->odd);
    *v_r = r1 + r2;
    *v_i = i1 + i2;
}
