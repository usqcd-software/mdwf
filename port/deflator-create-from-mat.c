#include <assert.h>
#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

/* FIXME total space occupied by umax will be NCV>NEV from Lanczos;
    1) realloc space 
        - what about alignment?
    2) reuse space for vmax vectors for EigCG 
        - will we realistically use EigCG after Lanczos?
*/

int
Q(create_deflator_inplace)(
        struct Q(Deflator) **deflator_ptr,
        struct QX(HalfFermionMat) **hfm_ptr, 
        int hfm_nev,
        int eigcg_vmax, int eigcg_nev,
        double eigcg_eps, int eigcg_umax)
{
    DECLARE_STATE;
    struct Q(Deflator) *d;
    int dim, Ls;

    CHECK_ARG0(*hfm_ptr);

    if (NULL == hfm_ptr || NULL == *hfm_ptr)
        return q(set_error)(state, 0, 
                "create_deflator_inplace_half_fermion_matrix(): NULL matrix");
    if (deflator_ptr == NULL)
        return q(set_error)(state, 0, 
                "create_deflator_inplace_half_fermion_matrix(): NULL pointer");

    d = q(malloc)(state, sizeof (struct Q(Deflator)));
    if (NULL == d)
        return q(set_error)(state, 0, 
                "create_deflator_inplace_half_fermion_matrix(): not enough memory");

    /* check data types */
#if defined(HAVE_LAPACK)
    {
        doublecomplex dc;
        assert( &(dc.r) == (double *)(&dc) &&
                &(dc.i) == (double *)(&dc) + 1 &&
                sizeof(dc) == 2 * sizeof(double));
    }
#elif defined(HAVE_GSL)
    {
        doublecomplex dc;
        assert( &(dc.r) == (double *)(&dc) &&
                &(dc.i) == (double *)(&dc) + 1 &&
                sizeof(dc) == 2 * sizeof(double));
        gsl_complex *gc = (gsl_complex *)(&dc);
        assert( &(dc.r) == &(gc->dat[0]) &&
                &(dc.i) == &(gc->dat[1]) &&
                sizeof(dc) == sizeof(*gc));
    }
#else 
#  error "no linear algebra library"
#endif

    /* TODO copy (*hfm_ptr)->m to d->u */

    /* TODO set EigCG to "full" condition */

    /* TODO set hfm_ptr->m to "empty" so that it is not destroyed */
    return 0;
}
