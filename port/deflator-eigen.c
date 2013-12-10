#include <mdwf.h>

#if defined(HAVE_LAPACK)
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  include <gsl/gsl_linalg.h>
#  include <gsl/gsl_sort_double.h>
#  include <gsl/gsl_sort_vector_double.h>
#else
#  error "no linear algebra library"
#endif

int
Q(deflator_eigen)(double *evals,
                  struct Q(Deflator) *d)
{
    int i;

    if ((d == 0) || (evals == 0))
        return 1;

    if (d->loading)
        return q(set_error)(d->state, 0, "deflator_eigen(): in loading state");
        
    for (i = 0; i < d->usize; i++)
        evals[i] = d->hevals[i];

    /* normal return */
    return 0; 
}
