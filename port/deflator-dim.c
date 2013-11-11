#include <mdwf.h>

int
Q(deflator_current_dim)(struct Q(Deflator) *d)
{
    if (d == NULL)
        return 0;
    return d->usize;
}
