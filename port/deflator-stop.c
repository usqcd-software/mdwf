#include <mdwf.h>

void
Q(deflator_stop)(struct Q(Deflator) *d)
{
    if (NULL == d)
        return;
    d->frozen = 1;
}
