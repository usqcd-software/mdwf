#include <mdwf.h>
#include <string.h>

int
Q(deflator_start_load)(struct Q(Deflator) *d)
{
  if (NULL == d)
    return 1;
  if (d->loading == 1)
    return 1;
  d->loading = 1;

  return 0;
}
