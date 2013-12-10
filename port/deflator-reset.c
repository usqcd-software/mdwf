#include <mdwf.h>
#include <string.h>

void
Q(deflator_reset)(struct Q(Deflator) *d)
{
    if (NULL == d)
        return;
    if (d->loading != 0 || ! d->do_eigcg)
      return;

    struct q(DeflatorEigcg) *d_e = &(d->df_eigcg);
    d_e->vsize = 0;
    memset(d_e->T, 0, d_e->vmax * d_e->vmax * sizeof(d_e->T[0]));
}
