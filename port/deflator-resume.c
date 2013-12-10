#include <mdwf.h>

void
Q(deflator_resume)(struct Q(Deflator) *d)
{
    if (NULL == d)
      return;
    if (d->loading != 0 || ! d->do_eigcg)
      return;

    if (d->usize < d->umax)
        d->df_eigcg.frozen = 0;
}
