#include <mdwf.h>

void
Q(deflator_resume)(struct Q(Deflator) *d)
{
    if (NULL == d)
        return;
    if (d->usize < d->umax)
        d->frozen = 0;
}
