#include <mdwf.h>
#include <string.h>

int
Q(deflator_stop_load)(struct Q(Deflator) *d)
{
  if (NULL == d)
    return 1;
  if (d->loading == 0)
    return 1;
  d->loading = 0;

  q(df_rebuild)(d);

  return 0;
}
