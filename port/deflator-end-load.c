#include <mdwf.h>
#include <string.h>

void
Q(deflator_end_load)(struct Q(Deflator) *d)
{
    if (NULL == d)
        return;
    if (d->loading == 0)
      return;
    q(df_rebuild)(d);
    d->loading = 0;
}
