#include <mdwf.h>

void
Q(deflator_stop)(struct Q(Deflator) *d)
{
    if (NULL == d)
      return;
    if (d->loading != 0 || ! d->do_eigcg)
      return;
    
    d->df_eigcg.frozen = 1;
}
