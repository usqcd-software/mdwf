#ifndef DEFLATOR_H_bzjvAey3yv1jyafMycej
#define DEFLATOR_H_bzjvAey3yv1jyafMycej

#include "deflator-la.h"

#if defined(HAVE_LAPACK)
/* FIXME include LAPACK c-prototypes */
#elif defined(HAVE_GSL)
#  include <gsl/gsl_vector.h>
#  include <gsl/gsl_matrix.h>
#  include <gsl/gsl_eigen.h>
#else
#  error "no linear algebra library"
#endif

/* sub-struct of deflator, allocated only for EigCG */
struct q(DeflatorEigcg) {
    int                 vmax;
    int                 vsize;
    int                 nev;
    int                 frozen;

    double              eps;
    double              resid_norm_sq_min;
    
    latmat_c            V;
    latmat_c            tmp_V;
    /* FIXME may be reduced to 1 aux vector */
    latvec_c            work_c_1;
    latvec_c            work_c_2;
    doublecomplex       *T;
    double              *hevals;
    doublecomplex       *hevecs2;
    long int            lwork;
    doublecomplex       *zwork;
#if defined(HAVE_LAPACK)
    doublecomplex       *hevecs1;
    doublecomplex       *tau;
    double              *rwork;
#elif defined(HAVE_GSL)
    gsl_matrix_complex  *gsl_T_full;
    gsl_matrix_complex  *gsl_hevecs1;
    gsl_vector          *gsl_hevals1;
    gsl_eigen_hermv_workspace *gsl_wkspace1;
    gsl_matrix_complex  *gsl_T_m1;              /* can reuse gsl_T_full */
    gsl_matrix_complex  *gsl_hevecs2;           /* can reuse gsl_hevecs1 */
    gsl_vector          *gsl_hevals2;           /* can reuse gsl_hevals1 */
    gsl_eigen_hermv_workspace *gsl_wkspace2;    /* can reuse gsl_wkspace1 */
    gsl_matrix_complex  *gsl_T_proj;            /* can reuse gsl_T_full */
    gsl_matrix_complex  *gsl_hevecs3;           /* can reuse gsl_T_hevecs1 */
    gsl_eigen_hermv_workspace *gsl_wkspace3;    /* can reuse gsl_T_wkspace1 */
    gsl_matrix_complex  *gsl_QR;                
    gsl_matrix_complex  *gsl_Q_unpack;
    gsl_matrix_complex  *gsl_tmp_MxS;           
    gsl_vector_complex  *gsl_tau;
    size_t              *hevals_select1;
    size_t              *hevals_select2;        /* can reuse hevals_select1 */
#else
#  error "no linear algebra library"
#endif
};

struct Q(Deflator) {
    struct Q(State) *state;

    /* non-zero iff the deflator has allocated EigCG workspace */
    int                 do_eigcg;   
    struct q(DeflatorEigcg) df_eigcg;

    int                 umax;       /* max size of U-space */
    int                 usize;      /* cur size of U-space */
    int                 loading;

    /* incr_eig current state */
    latmat_c            U;          /* deflator space */
    /* aux workspace lat.vectors */
    latvec_c            work_c_1;
    latvec_c            work_c_2;
    doublecomplex       *zwork;
    doublecomplex       *H;         /* op in U-space: U^dag.A.U */
    doublecomplex       *H_ev;      /* wkspace for e.val calculation */ 
    double              *hevals;    /* current U^dag.A.U eigenvalues */
    doublecomplex       *C;         /* Cholesky decomposition of H */
    
    /* aux workspace for eigenvalue calcs */
#if defined(HAVE_LAPACK)
    long int            eig_lwork;
    doublecomplex       *eig_zwork;
    double              *eig_rwork;
#elif defined(HAVE_GSL)
    gsl_eigen_herm_workspace *gsl_eig_wkspace;
#else
#  error "no linear algebra library"
#endif

};


int q(init_df_eigcg) (
        struct q(DeflatorEigcg) *df_eigcg, struct Q(State) *s, 
        int vmax, int nev, double eps);
void q(fini_df_eigcg) (struct q(DeflatorEigcg) *df_eigcg, struct Q(State) *s);

int q(init_deflator)(
        struct Q(Deflator) *df, struct Q(State) *s, 
        int umax, int do_eigcg, int vmax, int nev, double eps);
void q(fini_deflator)(struct Q(Deflator) *df, struct Q(State) *s);

int q(deflator_calc_evals)(struct Q(Deflator) *d);


#endif/*DEFLATOR_H_bzjvAey3yv1jyafMycej*/
