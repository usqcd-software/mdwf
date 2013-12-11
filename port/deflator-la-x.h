//#ifndef DEFLATOR_LA_X_H_Qplwonz8vdfqchXt7fyi
//#define DEFLATOR_LA_X_H_Qplwonz8vdfqchXt7fyi

#if defined(HAVE_LAPACK)
#  include <f2c_types.h>
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  undef  CHECK_GSL_STATUS
#  define CHECK_GSL_STATUS(func) do { int status__ = func; \
    if (status__) fprintf(stderr,"%s:%d: %s\n", __FILE__, __LINE__, gsl_strerror(status__)); \
    assert(0 == status__); } while(0)
#  include <gsl/gsl_vector.h>
#  include <gsl/gsl_matrix.h>
#  include <gsl/gsl_eigen.h>
#  ifndef DEFLATOR_H_HAVE_DCOMPLEX
#  define DEFLATOR_H_HAVE_DCOMPLEX
        typedef struct {
            double r, i;
        } doublecomplex;
#   endif/*DEFLATOR_H_HAVE_DCOMPLEX*/
#else
#  error "no linear algebra library"
#endif 


/***************************/
/*  deflator vector defs   */
/***************************/

struct Fermion;
typedef struct {
    int dim; /* aka size */
    int Ls;
    struct Fermion *f;
    void *mem_ptr;
    size_t mem_size;
} qx(defl_vec);

#undef  defl_vec_is_null
#define defl_vec_is_null(p) (NULL == (p)->f)
#undef  defl_vec_set_null
#define defl_vec_set_null(p) do { (p)->f = NULL; } while(0)

qx(defl_vec) qx(defl_vec_alloc)(struct Q(State) *state);
qx(defl_vec) qx(defl_vec_view)(struct Q(State) *state, struct Fermion *f);
void qx(defl_vec_free)(struct Q(State) *state, qx(defl_vec) *v);
void qx(defl_vec_zero)(qx(defl_vec) x); /* x <- 0 */
void qx(defl_vec_copy)(qx(defl_vec) x, qx(defl_vec) y); /* y <- x */
/* return flop count */
int qx(defl_vec_dotu)(doublecomplex *res, qx(defl_vec) x, qx(defl_vec) y);
int qx(defl_vec_scal)(double alpha, qx(defl_vec) x);
int qx(defl_vec_axpy)(double alpha, qx(defl_vec) x, qx(defl_vec) y);
int qx(defl_vec_nrm2)(double *res, qx(defl_vec) x);

/* lin. operator for deflator: y <- MxM(x) 
   flop counding is done through ws */
void qx(defl_vec_linop)(
        qx(defl_vec) y, 
        qx(defl_vec) x, 
        struct qx(MxM_workspace) *ws); 

/*******************************/
/*  deflator matrix defs       */
/*******************************/
struct vFermion;
typedef struct {
    int dim; /* aka size */
    int Ls;
    ptrdiff_t stride;
    int size;
    int begin;
    int len;
    struct vFermion *fv;
    void *mem_ptr;
    size_t mem_size;
} qx(defl_mat);

#undef  defl_mat_is_null
#define defl_mat_is_null(p) (NULL == (p)->fv)
#undef  defl_mat_set_null
#define defl_mat_set_null(p) do { (p)->fv = NULL; } while(0)

qx(defl_mat) qx(defl_mat_alloc)(struct Q(State) *state, int ncol);
qx(defl_mat) qx(defl_mat_view)(struct Q(State) *state, int size, struct vFermion *fv);
void qx(defl_mat_free)(struct Q(State) *state, qx(defl_mat) *m);
//int qx(defl_mat_convert_to_blas)(qx(defl_mat) *m);
//int qx(defl_mat_convert_from_blas)(qx(defl_mat) *m);
void qx(defl_mat_copy)(qx(defl_mat) m1, qx(defl_mat) m2); /* m2 <- m1 */
/* create a submatrix of subset of columns
   only a 'view' is created, no allocation is performed
   do not try to 'free' submatrix: it may result in memory error */
qx(defl_mat) qx(defl_mat_submat_col)(qx(defl_mat) m, int col, int ncol);
int qx(defl_mat_insert_col)(qx(defl_mat) m, int col, qx(defl_vec) v);
int qx(defl_mat_get_col)(qx(defl_mat) m, int col, qx(defl_vec) v);


/* C <- A^\dag * B, A:lat*m, B:lat*n, C:m*n */
int qx(defl_lmH_dot_lm)(int m, int n,
                        qx(defl_mat) a, 
                        qx(defl_mat) b, 
                        doublecomplex *c, int ldc);
/* y <- A^\dag x, A:lat*m, x:lat, y:m */
int qx(defl_lmH_dot_lv)(int m,
                        qx(defl_mat) a, 
                        qx(defl_vec) x, 
                        doublecomplex *y);
/* C <- A * B, A:lat*k, B:k*n, C:lat*n */
int qx(defl_lm_dot_zm)(int n, int k,
                       qx(defl_mat) a,
                       doublecomplex *b, int ldb, 
                       qx(defl_mat) c);
/* y <- A * x, A:lat*n, x:n, y:lat */
int qx(defl_lm_dot_zv)(int n,
                       qx(defl_mat) a, 
                       doublecomplex *x,
                       qx(defl_vec) y);


//#endif/*DEFLATOR_LA_X_H_Qplwonz8vdfqchXt7fyi*/
