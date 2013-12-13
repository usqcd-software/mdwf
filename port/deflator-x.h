//#ifndef DEFLATOR_H_bzjvAey3yv1jyafMycej
//#define DEFLATOR_H_bzjvAey3yv1jyafMycej


#include "deflator-la-x.h"



/* sub-struct of deflator, allocated only for EigCG */
struct qx(DeflatorEigcg) {
    int                 vmax;
    int                 vsize;
    int                 nev;
    int                 frozen;

    double              eps;
    double              resid_norm_sq_min;
    
    qx(defl_mat)        V;
    qx(defl_mat)        tmp_V;
    qx(defl_vec)        work_c_1;
    qx(defl_vec)        work_c_2;
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

struct QX(Deflator) {
    struct Q(State) *state;

    /* non-zero iff the deflator has allocated EigCG workspace */
    int                 do_eigcg;   
    struct qx(DeflatorEigcg) df_eigcg;

    int                 umax;       /* max size of U-space */
    int                 usize;      /* cur size of U-space */
    int                 loading;

    /* incr_eig current state */
    qx(defl_mat)        U;          /* deflator space */
    /* aux workspace lat.vectors */
    qx(defl_vec)        work_c_1;
    qx(defl_vec)        work_c_2;
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

void qx(init_defl_eigcg_set_null)(struct qx(DeflatorEigcg) *d_e);
int qx(init_defl_eigcg) (
        struct qx(DeflatorEigcg) *df_eigcg, 
        struct Q(State) *s, 
        int vmax, int nev, double eps);
int qx(init_deflator)(
        struct QX(Deflator) *df, struct Q(State) *s, 
        int umax, qx(defl_mat) *u, int usize,
        int do_eigcg, int vmax, int nev, double eps);


void qx(fini_defl_eigcg)(
        struct qx(DeflatorEigcg) *df_eigcg, 
        struct Q(State) *s);
void qx(fini_deflator)(
        struct QX(Deflator) *df, 
        struct Q(State) *s);


int qx(defl_ortho_uspace)(qx(defl_vec) vec,
                          struct QX(Deflator) *df,
                          int u_lo, int u_hi);
int qx(defl_solve_in_uspace)(struct QX(Deflator)    *df,
                             struct Fermion         *x,
                             struct Fermion         *b);
int qx(defl_calc_evals)(struct QX(Deflator) *d);
int qx(defl_inject)(struct QX(Deflator) *deflator,
                    struct qx(MxM_workspace) *ws,
                    int u_lo, int u_hi, int u_pos,
                    qx(defl_vec) vec);
int qx(defl_inject_back)(struct QX(Deflator)  *deflator,
                         struct qx(MxM_workspace) *ws,
                         qx(defl_vec) vec);
int qx(defl_recalc_mat)(struct QX(Deflator) *deflator,
                        struct qx(MxM_workspace) *ws);
int qx(defl_rebuild)(struct QX(Deflator) *deflator);


int qx(defl_eigcg_preamble)(struct QX(Deflator)       *deflator,
                            struct Fermion            *psi_e,
                            struct Fermion            *rho_e,
                            double                    *rho_norm2,
                            struct Fermion            *chi_e, /* const ! */
                            struct qx(MxM_workspace)  *ws,
                            unsigned int               options);
int qx(defl_eigcg_update0)(struct QX(Deflator)      *deflator,
                           double                    a1,
                           double                    b1,
                           double                    a0,
                           double                    b0,
                           double                    r,
                           struct Fermion           *rho,
                           unsigned int              options);
int qx(defl_eigcg_update1)(struct QX(Deflator)      *deflator,
                           double                    a1,
                           double                    b1,
                           double                    a0,
                           double                    b0,
                           double                    r,
                           struct Fermion           *rho,
                           struct Fermion           *A_rho,
                           unsigned int              options);
int qx(defl_eigcg_postamble)(struct QX(Deflator)       *deflator,
                             struct qx(MxM_workspace)  *ws,
                             unsigned int               options);

//#endif/*DEFLATOR_H_bzjvAey3yv1jyafMycej*/
