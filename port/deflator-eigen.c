#include <mdwf.h>

int
Q(deflator_eigen)(double *eigen_values,
                  struct Q(Deflator) *d)
{
    if ((d == 0) || (eigen_values == 0))
        return 1;
    
    /* XXX export deflator eigenvalues */
    return 1; /* no eigenvalues so far */
}
