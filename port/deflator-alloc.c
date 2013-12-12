#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>
#include <assert.h>

#define ds      sizeof(double)
#define zs      sizeof(doublecomplex)


/* fill & allocate values in `d_e' (`d_e' must be allocated) */
int 
qx(init_defl_eigcg) (
        struct qx(DeflatorEigcg) *d_e, 
        struct Q(State) *s,
        int vmax, int nev, double eps)
{
    if (s == NULL || s->error_latched)
        return 1;

    if (d_e == NULL)
        return q(set_error)(s, 0, "init_df_eigcg(): NULL pointer");
  

    /* first, set all to null */
    defl_mat_set_null(&(d_e->V));
    defl_mat_set_null(&(d_e->tmp_V));
    defl_vec_set_null(&(d_e->work_c_1));
    defl_vec_set_null(&(d_e->work_c_2));
    d_e->T                = NULL;
    d_e->lwork            = 0;
    d_e->zwork            = NULL;
    d_e->hevals           = NULL;
    d_e->hevecs2          = NULL;
#if defined(HAVE_LAPACK)
    d_e->hevecs1          = NULL;
    d_e->tau              = NULL;
    d_e->rwork            = NULL;
#elif defined(HAVE_GSL)
    d_e->gsl_T_full       = NULL;
    d_e->gsl_hevecs1      = NULL;
    d_e->gsl_hevals1      = NULL;
    d_e->gsl_wkspace1     = NULL;
    d_e->gsl_T_m1         = NULL;
    d_e->gsl_hevecs2      = NULL;
    d_e->gsl_hevals2      = NULL;
    d_e->gsl_wkspace2     = NULL;
    d_e->gsl_T_proj       = NULL;
    d_e->gsl_hevecs3      = NULL;
    d_e->gsl_wkspace3     = NULL;
    d_e->gsl_QR           = NULL;
    d_e->gsl_Q_unpack     = NULL;
    d_e->gsl_tmp_MxS      = NULL;
    d_e->gsl_tau          = NULL;
    d_e->hevals_select1   = NULL;
    d_e->hevals_select2   = NULL;
#else
#  error "no linear algebra library"
#endif

    /* allocate */
#define ds      sizeof(double)
#define zs      sizeof(doublecomplex)
    d_e->V                = qx(defl_mat_alloc)(s, vmax);
    d_e->tmp_V            = qx(defl_mat_alloc)(s, 2*nev);
    d_e->work_c_1         = qx(defl_vec_alloc)(s);
    d_e->work_c_2         = qx(defl_vec_alloc)(s);
    d_e->T                = q(malloc)(s, vmax * vmax * zs);
    d_e->hevals           = q(malloc)(s, vmax * ds);
    d_e->hevecs2          = q(malloc)(s, vmax * vmax * zs);
    d_e->lwork            = 2*vmax;
    d_e->zwork            = q(malloc)(s, d_e->lwork * zs);
#if defined(HAVE_LAPACK)
    d_e->hevecs1          = q(malloc)(s, vmax * vmax * zs);
    d_e->tau              = q(malloc)(s, vmax * zs);
    d_e->rwork            = q(malloc)(s, 3 * vmax * ds);
#elif defined(HAVE_GSL)
    d_e->gsl_T_full       = gsl_matrix_complex_alloc(vmax, vmax);
    d_e->gsl_hevecs1      = gsl_matrix_complex_alloc(vmax, vmax);
    d_e->gsl_hevals1      = gsl_vector_alloc(vmax);
    d_e->gsl_wkspace1     = gsl_eigen_hermv_alloc(vmax);
    d_e->gsl_T_m1         = gsl_matrix_complex_alloc(vmax-1, vmax-1);
    d_e->gsl_hevecs2      = gsl_matrix_complex_alloc(vmax-1, vmax-1);
    d_e->gsl_hevals2      = gsl_vector_alloc(vmax-1);
    d_e->gsl_wkspace2     = gsl_eigen_hermv_alloc(vmax-1);
    d_e->gsl_T_proj       = gsl_matrix_complex_alloc(2*nev, 2*nev);
    d_e->gsl_hevecs3      = gsl_matrix_complex_alloc(2*nev, 2*nev);
    d_e->gsl_wkspace3     = gsl_eigen_hermv_alloc(2*nev);
    d_e->gsl_QR           = gsl_matrix_complex_alloc(vmax, 2*nev);
    d_e->gsl_Q_unpack     = gsl_matrix_complex_alloc(vmax, vmax);
    d_e->gsl_tmp_MxS      = gsl_matrix_complex_alloc(vmax, 2*nev);
    d_e->gsl_tau          = gsl_vector_complex_alloc(2*nev);
    d_e->hevals_select1   = q(malloc)(s, vmax * sizeof(d_e->hevals_select1[0]));
    d_e->hevals_select2   = q(malloc)(s, vmax * sizeof(d_e->hevals_select2[0]));
#else
#  error "no linear algebra library"
#endif

    /* check allocation */
    if (       defl_mat_is_null(&(d_e->V))
            || defl_mat_is_null(&(d_e->tmp_V))
            || defl_vec_is_null(&(d_e->work_c_1))
            || defl_vec_is_null(&(d_e->work_c_2))
            || NULL == d_e->T
            || NULL == d_e->hevals
            || NULL == d_e->hevecs2
            || NULL == d_e->zwork
#if defined(HAVE_LAPACK)
            || NULL == d_e->hevecs1
            || NULL == d_e->tau
            || NULL == d_e->rwork
#elif defined(HAVE_GSL)
            || NULL == d_e->gsl_T_full
            || NULL == d_e->gsl_hevecs1
            || NULL == d_e->gsl_hevals1
            || NULL == d_e->gsl_wkspace1
            || NULL == d_e->gsl_T_m1
            || NULL == d_e->gsl_hevecs2
            || NULL == d_e->gsl_hevals2
            || NULL == d_e->gsl_wkspace2
            || NULL == d_e->gsl_T_proj
            || NULL == d_e->gsl_hevecs3
            || NULL == d_e->gsl_wkspace3
            || NULL == d_e->gsl_QR
            || NULL == d_e->gsl_Q_unpack
            || NULL == d_e->gsl_tmp_MxS
            || NULL == d_e->gsl_tau
            || NULL == d_e->hevals_select1
            || NULL == d_e->hevals_select2
#else
#  error "no linear algebra library"
#endif
            ) {
        qx(fini_defl_eigcg)(d_e, s);
        return q(set_error)(s, 0, "init_df_eigcg(): not enough memory");
    }

    BEGIN_TIMING(s);

    d_e->vmax   = vmax;
    d_e->vsize  = 0;
    d_e->nev    = nev;
    d_e->eps    = eps;
    d_e->frozen = 0;

    END_TIMING(s, 0, 0, 0);
    /* TODO check that everything is allocated; otherwise, call 'free' */

  return 0;
}


/* fill & allocate values in `df' (`df' must be allocated) ;
   latmat `u' is set to null on return to prevent deallocation or gc */
int
qx(init_deflator)(
        struct QX(Deflator) *df, 
        struct Q(State) *s,
        int umax, qx(defl_mat) *u, int usize,
        int do_eigcg, int vmax, int nev, double eps)
{
    int status = 0;
    if (s == NULL || s->error_latched)
        return 1;

    if (df == NULL)
        return q(set_error)(s, 0, "init_deflator(): NULL pointer");
    
    df->do_eigcg = do_eigcg;
    if (do_eigcg) {
        if (0 != (status = qx(init_defl_eigcg)(&(df->df_eigcg), s, vmax, nev, eps)))
            return status;
    }
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


    /* set all to NULL */
    defl_mat_set_null(&(df->U));
    defl_vec_set_null(&(df->work_c_1));
    defl_vec_set_null(&(df->work_c_2));
    df->zwork       = NULL;
    df->H           = NULL;
    df->H_ev        = NULL;
    df->hevals      = NULL;
    df->C           = NULL;
#if defined(HAVE_LAPACK)
    df->eig_lwork   = 0;
    df->eig_zwork   = NULL;
    df->eig_rwork   = NULL;
#elif defined(HAVE_GSL)
    df->gsl_eig_wkspace = NULL;
#else
#  error "no linear algebra library"
#endif
    
    /* allocate (or load external) */
    if (NULL == u) {
        df->U       = qx(defl_mat_alloc)(s, umax);
        df->umax    = umax;
        df->usize   = 0;
    } else {
        int u_len = u->len;
        if (u_len < usize) {
            qx(fini_deflator)(df, s);
            return q(set_error)(s, 0, "init_deflator(): matrix is smaller than USIZE");
        }
        if (u_len < umax) {
            qx(fini_deflator)(df, s);
            return q(set_error)(s, 0, "init_deflator(): UMAX is larger than provided matrix");
        }
        if (umax < usize)
            umax = usize;

        df->umax    = umax;
        df->usize   = usize;
        qx(defl_mat_swap)(u, &(df->U));
    }

    df->work_c_1    = qx(defl_vec_alloc)(s);
    df->work_c_2    = qx(defl_vec_alloc)(s);
    df->zwork       = q(malloc)(s, umax * zs);
    df->H           = q(malloc)(s, umax * umax * zs);
    df->H_ev        = q(malloc)(s, umax * umax * zs);
    df->hevals      = q(malloc)(s, umax * ds);
    df->C           = q(malloc)(s, umax * umax * zs);
#if defined(HAVE_LAPACK)
    df->eig_lwork   = 2 * umax;
    df->eig_zwork   = q(malloc)(s, df->lwork * zs);
    df->eig_rwork   = q(malloc)(s, 3 * umax * ds);
#elif defined(HAVE_GSL)
    df->gsl_eig_wkspace = gsl_eigen_herm_alloc(umax);
#else
#  error "no linear algebra library"
#endif

    /* check allocation */
    if (       defl_mat_is_null(&(df->U))
            || defl_vec_is_null(&(df->work_c_1))
            || defl_vec_is_null(&(df->work_c_2))
            || NULL == df->zwork
            || NULL == df->H
            || NULL == df->H_ev
            || NULL == df->hevals
            || NULL == df->C
#if defined(HAVE_LAPACK)
            || NULL == df->eig_zwork
            || NULL == df->eig_rwork
#elif defined(HAVE_GSL)
            || NULL == df->gsl_eig_wkspace
#else
#  error "no linear algebra library"
#endif
            ) {
        qx(fini_deflator)(df, s);
        return q(set_error)(s, 0, "init_deflator(): not enough memory");
    }

    df->state   = s;
    df->loading = 0;
    return 0;
}


int
QX(create_deflator)(struct QX(Deflator) **deflator_ptr,
                   struct Q(State) *s,
                   int vmax, int nev, 
                   double eps, int umax)
{
    struct QX(Deflator) *d;
    int status;
    if (s == NULL || s->error_latched)
        return 1;

    if (deflator_ptr == NULL)
        return q(set_error)(s, 0, "create_deflator(): NULL pointer");
  
    *deflator_ptr = NULL;
    d = q(malloc)(s, sizeof (struct QX(Deflator)));
    if (d == 0)
        return q(set_error)(s, 0, "allocate_deflator(): not enough memory");

    BEGIN_TIMING(s);
    if (0 != (status = qx(init_deflator)(d, s, 
                                         umax, NULL, 0,
                                         1/*do_eigcg*/, vmax, nev, eps)))
        return status;
    END_TIMING(s, 0, 0, 0);

    *deflator_ptr = d;

    return 0;
}
