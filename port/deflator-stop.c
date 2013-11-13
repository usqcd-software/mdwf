#include <mdwf.h>

void
Q(deflator_stop)(struct Q(Deflator) *d)
{
    if (NULL == d)
      return;
    if (d->loading != 0)
      return;
    d->frozen = 1;
}
