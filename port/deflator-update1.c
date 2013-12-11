#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>
#include <math.h>
#include <qmp.h>

#if defined(HAVE_LAPACK)
#  include <blas.h>
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  include <gsl/gsl_linalg.h>
#  include <gsl/gsl_eigen.h>
#  include <gsl/gsl_blas.h>
#  include <gsl/gsl_complex_math.h>
#  include <gsl/gsl_sort_double.h>
#  include <gsl/gsl_sort_vector_double.h>
#else
#  error "no linear algebra library"
#endif

/* XXX uses WK: df->d_e.work_c_1 and many other vars in df->d_e 
        reserved for update1 */
int
qx(defl_update1)(
        struct QX(Deflator)     *df,
        double                   alpha, 
        double                   beta, 
        double                   alpha_prev, 
        double                   beta_prev, 
        double                   resid_norm_sq, 
        struct Fermion          *resid,
        struct Fermion          *A_resid,
        unsigned int             options)
{
    int i, j;
    long int vmax, vsize, nev;
    struct qx(DeflatorEigcg) *d_e = NULL;
#if defined(HAVE_LAPACK)
    char cV = 'V',
         cU = 'U',
         cL = 'L',
         cR = 'R',
         cC = 'C',
         cN = 'N';
    long int info = 0;
    long int tmp_i;
#endif

    if (NULL == df || 
            NULL == df->state ||
            df->df_eigcg.frozen)
        return 0;       /* skip call */
    d_e = &(df->df_eigcg);

    /* FIXME allow updates with vsize < vmax; 
       what about the final CG iteration? are (vsize-2*nev) vectors discarded? */
    if (d_e->vsize < d_e->vmax) 
        return -1;      /* should call update0 instead */
    assert(d_e->vmax == d_e->vsize && 1 < d_e->vsize);

    /* set [vsize-1, vsize-1] elem */
    doublecomplex *pT = d_e->T + (d_e->vsize - 1 ) * (1 + d_e->vmax);
    pT->r = 1. / alpha + beta_prev / alpha_prev;
    pT->i = 0.0;

    /* do restart */
    nev         = d_e->nev;
    d_e->vsize  = 2 * nev;
    vmax        = d_e->vmax;
    vsize       = d_e->vsize;

#if defined(HAVE_LAPACK)
    char cV = 'V',
         cU = 'U',
         cL = 'L',
         cR = 'R',
         cC = 'C',
         cN = 'N';
    long int info = 0;
    long int tmp_i;


    /* diagonalize T:vmax*vmax matrix 
       requires zwork size = lwork >= 2*vmax-1, rwork size >=3*vmax-2 */
    /* w,X <- eigh(T)  # w:=hevals, X:=hevecs1, T:=hevecs1 */
    memcpy(d_e->hevecs1, d_e->T, vmax * vmax * sizeof(d_e->T[0]));
    zheev_(&cV, &cU, &vmax, d_e->hevecs1, &vmax, 
           d_e->hevals, d_e->zwork, &(d_e->lwork), 
           d_e->rwork, &info, 1, 1);
    assert(0 == info);
    if (options & QOP_MDWF_LOG_EIG_UPDATE1) {
        int i;
        for (i = 0; i < vmax; i++)
            qf(zprint)(df->state, "update1", "T0 %4d %17.9e", i, d_e->hevals[i]);
                       
    }

    /* diagonalize T:(vmax-1)*(vmax-1) matrix
       requires zwork size = lwork >= 2*vmax-1 */
    /* w',X' <- eigh(T[:-1,:-1] # w':=hevals, X':=hevecs2, T:=hevecs2*/
    memcpy(d_e->hevecs2, d_e->T, vmax * vmax * sizeof(d_e->T[0]));
    tmp_i = vmax - 1;
    zheev_(&cV, &cU, &tmp_i, d_e->hevecs2, &vmax,
           d_e->hevals, d_e->zwork, &(d_e->lwork),
           d_e->rwork, &info, 1, 1);
    assert(0 == info);

    /* select first nev vectors from both spaces */
    /* Q'[:,:2*nev] <- concat_in_dim2[ X[:,:nev], X'[:,:nev] ] # Q':=hevecs1 */
    memcpy(d_e->hevecs1 + d_e->nev * vmax, d_e->hevecs2, 
           d_e->nev * vmax * sizeof(d_e->hevecs1[0]));
    /* fill [vmax-1, nev:2*nev] with zeros */
    for (i = d_e->nev; i < 2 * d_e->nev; i++) {
        doublecomplex *p = d_e->hevecs1 + (i + 1) * vmax - 1;
        p->r = 0.0;
        p->i = 0.0;
    }

    /* QR factorization/orthogonalization of 2*nev columns
       requires zwork size = lwork >= vsize=2*nev(?), tau size >= 2*nev */
    /* Q,R <- qr_decomp(Q') # Q:=L{hevecs1},tau, R:=U{hevecs1} */
    zgeqrf_(&vmax, &vsize, d_e->hevecs1, &vmax, 
            d_e->tau, d_e->zwork, &(d_e->lwork),
            &info);
    assert(0 == info);

    /* compute hevecs2 <- Q^H T Q */
    memcpy(d_e->hevecs2, d_e->T, vmax * vmax * sizeof(d_e->T[0]));
    /* requires size of zwork = lwork >= vmax ; optimal-?*/
    /* Tp <- Q^H.T.Q  # Q:=L{hevecs1},tau, T:=hevecs2 */
    zunmqr_(&cL, &cC, &vmax, &vmax, &vsize, 
            d_e->hevecs1, &vmax, d_e->tau,
            d_e->hevecs2, &vmax, d_e->zwork,
            &(d_e->lwork), &info, 1, 1);
    assert(0 == info);
    zunmqr_(&cR, &cN, &vsize, &vmax, &vsize, 
            d_e->hevecs1, &vmax, d_e->tau,
            d_e->hevecs2, &vmax, d_e->zwork,
            &(d_e->lwork), &info, 1, 1);
    assert(0 == info);

    /* compute eigenpairs of Q^H T Q = Z M Z^H, hevecs2 <- Z */
    /* w'',Z <- eigh(Tp) # w'':=hevals, Z:=hevecs2, Tp:=hevecs2 */
    zheev_(&cV, &cU, &vsize, d_e->hevecs2, &vmax, 
           d_e->hevals, d_e->zwork, &(d_e->lwork),
           d_e->rwork, &info, 1, 1);
    assert(0 == info);

    /* fill Z[2nev:vmax, 0:2nev] with zeros */
    for (j = 0; j < vsize; j++) {
        doublecomplex *p = d_e->hevecs2 + j * vmax;
        for (i = vsize ; i < vmax; i++)
            p[i].r = p[i].i = 0.0;
    }
    /* compute hevecs2 <- Q Z */
    /* (QZ) <- Q.Z # Q:=L{hevecs1},tau, Z:=hevecs2 */
    zunmqr_(&cL, &cN, &vmax, &vsize, &vsize,
            d_e->hevecs1, &vmax, d_e->tau,
            d_e->hevecs2, &vmax, d_e->zwork,
            &(d_e->lwork), &info, 1, 1);
    assert(0 == info);
    

#elif defined(HAVE_GSL)

    gsl_matrix_complex_view gsl_T = gsl_matrix_complex_view_array(
            (double *)d_e->T, 
            vmax, vmax);
    gsl_matrix_complex_transpose(&gsl_T.matrix);    /* d_e->T uses FORTRAN matrix indexing */
    /* eigenpairs of T */
    /* gsl_hevals1, gsl_hevecs1 <- eigh(gsl_T) # gsl_T_full:=WK */
    gsl_matrix_complex_memcpy(d_e->gsl_T_full, &gsl_T.matrix);
    CHECK_GSL_STATUS(gsl_eigen_hermv(
            d_e->gsl_T_full, 
            d_e->gsl_hevals1, 
            d_e->gsl_hevecs1,
            d_e->gsl_wkspace1));
    CHECK_GSL_STATUS(gsl_sort_smallest_index(d_e->hevals_select1, d_e->nev, 
            gsl_vector_const_ptr(d_e->gsl_hevals1, 0), 1, vmax));

    if (options & QOP_MDWF_LOG_EIG_UPDATE1) {
        int i, j;
        
        for (i = 0; i < d_e->nev; i++) {
            j = d_e->hevals_select1[i];
            qf(zprint)(df->state, "update1", "T0 %4d %17.9e", i,
                       gsl_vector_get(d_e->gsl_hevals1, j));
        }
    }

    /* eigenpairs of T[:-1, :-1] */
    /* gsl_hevals2, gsl_hevecs2 <- eigh(gsl_T) # gsl_T_m1:=WK */
    gsl_matrix_complex_view gsl_T_m1 = gsl_matrix_complex_submatrix(
            &gsl_T.matrix, 
            0, 0, vmax-1, vmax-1);
    CHECK_GSL_STATUS(gsl_matrix_complex_memcpy(
                d_e->gsl_T_m1, 
                &gsl_T_m1.matrix));
    CHECK_GSL_STATUS(gsl_eigen_hermv(
            d_e->gsl_T_m1, 
            d_e->gsl_hevals2,
            d_e->gsl_hevecs2,
            d_e->gsl_wkspace2));
    CHECK_GSL_STATUS(gsl_sort_smallest_index(d_e->hevals_select2, d_e->nev,
            gsl_vector_const_ptr(d_e->gsl_hevals2, 0), 1, vmax-1));

    /* construct Q = (Y[:nev], Y1[:nev]) */
    /* gsl_QR <- concat_in_dim2[X[:nev], X'[:nev]] */
    for (j = 0; j < d_e->nev; j++) {
        int j1 = d_e->hevals_select1[j];
        for (i = 0; i < vmax; i++)
            gsl_matrix_complex_set(d_e->gsl_QR, i, j,
                    gsl_matrix_complex_get(d_e->gsl_hevecs1, i, j1));
    }
    for (j = 0; j < d_e->nev; j++) {
        int j2 = d_e->hevals_select2[j];
        for (i = 0; i < vmax-1; i++)
            gsl_matrix_complex_set(d_e->gsl_QR, i, d_e->nev + j,
                    gsl_matrix_complex_get(d_e->gsl_hevecs2, i, j2));
        gsl_matrix_complex_set(d_e->gsl_QR, vmax-1, d_e->nev + j,
                gsl_complex_rect(0., 0.));
    }

    /* QR decomp */
    /* gsl_Q, gsl_tmp_MxS <- Q,R = qr_decomp(gsl_QR) */
    /* FIXME rewrite this section and below using implicit QR decomp
       without constructing the full matrix Q */
    CHECK_GSL_STATUS(gsl_linalg_complex_QR_decomp(
            d_e->gsl_QR, d_e->gsl_tau));
    CHECK_GSL_STATUS(gsl_linalg_complex_QR_unpack(
            d_e->gsl_QR, d_e->gsl_tau, 
            d_e->gsl_Q_unpack, d_e->gsl_tmp_MxS));
    gsl_matrix_complex_view gsl_Q = gsl_matrix_complex_submatrix(
            d_e->gsl_Q_unpack,
            0, 0, vmax, vsize);

    /* projecting matrix T on Q: Q^H . T . Q */
    /* gsl_T_proj <- gsl_Q^H . gsl_T . gsl_Q */
    CHECK_GSL_STATUS(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
            GSL_COMPLEX_ONE, &gsl_T.matrix, &gsl_Q.matrix,
            GSL_COMPLEX_ZERO, d_e->gsl_tmp_MxS));
    CHECK_GSL_STATUS(gsl_blas_zgemm(CblasConjTrans, CblasNoTrans,
            GSL_COMPLEX_ONE, &gsl_Q.matrix, d_e->gsl_tmp_MxS,
            GSL_COMPLEX_ZERO, d_e->gsl_T_proj));

    /* eigenpairs of Q^H . T . Q */
    /* gevals, gsl_hevecs3 <- eigh(gsl_T_proj) # gsl_T_proj:WK */
    gsl_vector_view gsl_hevals = gsl_vector_view_array(
            d_e->hevals, vsize);
    CHECK_GSL_STATUS(gsl_eigen_hermv(
            d_e->gsl_T_proj,
            &gsl_hevals.vector,
            d_e->gsl_hevecs3,
            d_e->gsl_wkspace3));
    /**/CHECK_GSL_STATUS(gsl_eigen_hermv_sort(
                &gsl_hevals.vector, d_e->gsl_hevecs3, GSL_EIGEN_SORT_VAL_ASC));
    
    /* compute rotation and transpose (FORTRAN matrix storage conventions) */
    /* hevecs2[F] <- (gsl_Q . gsl_hevecs3[C])^T */
    CHECK_GSL_STATUS(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
            GSL_COMPLEX_ONE, &gsl_Q.matrix, d_e->gsl_hevecs3,
            GSL_COMPLEX_ZERO, d_e->gsl_tmp_MxS));

    gsl_matrix_complex_view gsl_QZ_transp = gsl_matrix_complex_view_array_with_tda(
            (double *)d_e->hevecs2, vsize, vmax, vmax);
    CHECK_GSL_STATUS(gsl_matrix_complex_transpose_memcpy( 
                &gsl_QZ_transp.matrix, d_e->gsl_tmp_MxS));


#else
#  error "no linear algebra library"
#endif
    
    /* XXX broadcast rotation matrix and eigenvalues from the master node
       the rotation matrix may differ by reflection: 
       V -> V\times diag{.. \pm 1 .. } */
    QMP_broadcast((void *)(d_e->hevecs2), vsize * vmax * sizeof(d_e->hevecs2[0]));
    /* FIXME eigenvalues should be the same, and they are not used later; omit? */
    QMP_broadcast((void *)(d_e->hevals), vsize * sizeof(d_e->hevals[0]));

    /* rotate V[:, 0:vmax] space with (Q Z) */
    /* FIXME change dot in V . (QZ) into rotation by Q and by Z 
       to avoid using tmp_V */
    qx(defl_mat) tmp_V = q(defl_mat_submat_col)(d_e->tmp_V, 0, vsize);
    q(defl_lm_dot_zm)(vsize, vmax, 
                  d_e->V,
                  d_e->hevecs2, vmax, 
                  tmp_V);
    qx(defl_mat_copy)(tmp_V, q(defl_mat_submat_col)(d_e->V, 0, vsize));
    
    /* check eig convergence */
    if (resid_norm_sq < d_e->resid_norm_sq_min)
        return 3;   /* eig converged */
    double resid_norm = sqrt(resid_norm_sq);
    
    /* compute new T */
    memset(d_e->T, 0, vmax * vmax * sizeof(d_e->T[0]));
    qx(defl_lmH_dot_lv)(vsize, 
                      tmp_V, 
                      qx(defl_vec_view)(df->state, A_resid),
                      d_e->T + vsize * vmax);
    for (i = 0 ; i < vsize ; i++) {
        /* T[i,i] <- hevals[i] + 0*I */
        d_e->T[i * (vmax + 1)].r      = d_e->hevals[i];
        d_e->T[i * (vmax + 1)].i      = 0.0;
        
        /* T[:,vsize] <- V^H . A . r */
        d_e->T[i + vsize * vmax].r    /= resid_norm;
        d_e->T[i + vsize * vmax].i    /= resid_norm;
        
        /* T[vsize,:] <- T[:,vsize]^* */
        d_e->T[vsize + i * vmax].r    =  d_e->T[i + vsize * vmax].r;
        d_e->T[vsize + i * vmax].i    = -d_e->T[i + vsize * vmax].i;
    }

    /* remember the vector ||resid|| */
    qx(defl_vec) cur_r = d_e->work_c_1;
    /* V[:,vsize] <- ||resid|| */
    qx(defl_vec_copy)(qx(defl_vec_view)(df->state, resid), 
                      cur_r);
    qx(defl_vec_scal)(1. / resid_norm, cur_r);
    qx(defl_mat_insert_col)(d_e->V, d_e->vsize, cur_r);

    d_e->vsize += 1;

    return 0; /* normal */
}
