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

    /* matrix is destroyed by the eigenvalue algorithm */
    /* FIXME only triangular part is destroyed, so no extra storage 
       is necessary; the full matrix can be restored after computing evals */
    memcpy(d->H_ev, d->H, d->usize * d->umax * sizeof(d->C[0]));
#if HAVE_LAPACK
    {
        long int usize  = d->usize,
                 umax   = d->umax,
                 info   = 0;
        char cU = 'U',
            cN = 'N';
        /* FIXME use debug workspace here; perhaps, rename the fields? */
        zheev_(&cN, &cU, &usize, d->H_ev, &umax, d->debug_hevals, 
               d->debug_zwork, &(d->debug_lwork), d->debug_rwork,
               &info, 1, 1);
        if (0 != info)
            goto clearerr_1;
        for (i = 0; i < d->usize; i++)
            evals[i] = d->debug_hevals[i];
    }
#elif HAVE_GSL
    {
        gsl_vector_view gsl_hevals = gsl_vector_subvector(
            d->debug_gsl_hevals, 0, d->usize);
        gsl_matrix_complex_view gsl_C =
            gsl_matrix_complex_view_array_with_tda(
                (double *)d->H_ev, d->usize, d->usize, d->umax);

        /* transpose to convert matrix Fortran[row][col] -> C[col][row] */
        if (gsl_matrix_complex_transpose(&gsl_C.matrix))
            goto clearerr_1;
        if (gsl_eigen_herm(&gsl_C.matrix, &gsl_hevals.vector,
                                        d->debug_gsl_wkspace))
            goto clearerr_1;
        gsl_sort_vector(&gsl_hevals.vector);
        for (i = 0; i < d->usize; i++)
            evals[i] = gsl_vector_get(&gsl_hevals.vector, i);
    }
#else
#  error "no linear algebra library"
    return 1;
#endif
    /* normal return */
    return 0; 

    /* cleanup & error exit*/
clearerr_1:
    return 1;
}
