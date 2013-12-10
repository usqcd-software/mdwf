#define QOP_MDWF_DEFAULT_PRECISION 'F'
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

int q(deflator_calc_evals)(struct Q(Deflator) *d) 
{
    long int usize  = d->usize,
             umax   = d->umax;
    /* matrix is destroyed by the eigenvalue algorithm; need to copy */
    /* FIXME only triangular part is destroyed, so no extra storage 
       is necessary; the full matrix can be restored after computing evals */
    memcpy(d->H_ev, d->H, d->usize * d->umax * sizeof(d->H[0]));
#if HAVE_LAPACK
    {
        int info   = 0;
        char cU = 'U',
             cN = 'N';
        zheev_(&cN, &cU, &usize, d->H_ev, &umax, d->hevals, 
               d->eig_zwork, &(d->eig_lwork), d->eig_rwork,
               &info, 1, 1);
        assert(0 == info);
    }
#elif HAVE_GSL
    {
        gsl_vector_view gsl_hevals = gsl_vector_view_array(d->hevals, usize);
        gsl_matrix_complex_view gsl_H_ev = gsl_matrix_complex_view_array_with_tda(
                (double *)d->H_ev, usize, usize, d->umax);
        /* transpose to convert matrix Fortran[row][col] -> C[col][row] */
        CHECK_GSL_STATUS(gsl_matrix_complex_transpose(&gsl_H_ev.matrix));

        CHECK_GSL_STATUS(gsl_eigen_herm(&gsl_H_ev.matrix, &gsl_hevals.vector,
                                        d->gsl_eig_wkspace));
        gsl_sort_vector(&gsl_hevals.vector);
    }
#else
#  error "no linear algebra library"
#endif
    return 0;
}
