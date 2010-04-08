#include <mdwf.h>
#include <string.h>

void
Q(deflator_reset)(struct Q(Deflator) *d)
{
    if (NULL == d)
        return;
    d->vsize = 0;
    memset(d->T, 0, d->vmax * d->vmax * sizeof(d->T[0]));
}
