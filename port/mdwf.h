# ifndef QOP_MDWF_DEFAULT_PRECISION
#  define QOP_MDWF_DEFAULT_PRECISION 'D'
# endif

#ifndef MARK_B9BA8123_0F1A_40FD_8827_42266FE32F3E
#define MARK_B9BA8123_0F1A_40FD_8827_42266FE32F3E

# include <qop-mdwf3.h>
# include <stdlib.h>
# include <string.h>
# include <qmp.h>
# include <sys/time.h>

# define q(x) qop_mdwf_##x
# define qf(x) qop_f3_mdwf_##x
# define qd(x) qop_d3_mdwf_##x
# define Q(x) QOP_MDWF_##x
# define QF(x) QOP_F3_MDWF_##x
# define QD(x) QOP_D3_MDWF_##x

/* Cache size */
#define CACHE_LINE_SIZE 128
#define ALIGN(p,d) ((void *)((((ptrdiff_t)(p))+(d)+CACHE_LINE_SIZE-1) &  \
                            ~(CACHE_LINE_SIZE-1)))

/* QCD types (qa0 controls these definitions) */
struct SUnD;
struct SUnF;
struct FermionF;
struct FermionD;
struct VectorFermionF;
struct VectorFermionD;
struct ProjectedFermionF;
struct ProjectedFermionD;
struct MxM_workspaceD;
struct MxM_workspaceF;

/* Internal types */
struct local {
  int lo[Q(DIM)];
  int hi[Q(DIM)];
  int dx[Q(DIM)];
};

/* structs neighbor and up_pack are defined by qa0 */
struct neighbor;
struct up_pack;
struct down_pack;

struct eo_lattice {
  struct Q(State) *state;                  /* back pointer to state */
  int              face_size;              /* 4-d size of the face */
  int              body_size;              /* 4-d size of the body */
  int              full_size;              /* face + body */
  int             *lx2v;                   /* 4-d layout 2 vector translation */
  int             *v2lx;                   /* only for init */
  int              Ls;                     /* Ls */
  struct local    *local;                  /* points to state.local */

  struct neighbor *neighbor;               /* neighbor data (body,face) */
  int              send_up_size[Q(DIM)];   /* 4-d send size in each up-dir */
  struct up_pack  *up_pack[Q(DIM)];        /* 4-d (U,f) for up-face packing */
  int              send_down_size[Q(DIM)]; /* 4-d send size in each down-dir */
  struct down_pack *down_pack[Q(DIM)];      /* 4-d (f) for down-face packing */
  int              receive_up_size[Q(DIM)]; /* 4-d (U,f) up receive size */
  int              receive_down_size[Q(DIM)]; /* 4-d (f) down receive size */

  int              real_size;              /* 0, 4 or 8 */ 
  int              h_valid;                /* is .handle valid? */
  QMP_msghandle_t  handle;                 /* global send&receive handle */
  QMP_msghandle_t  th[4*Q(DIM)];           /* transitody handles */
  int              th_count;               /* number of valid th[] */
  QMP_msgmem_t     mh[4*Q(DIM)];           /* memory handles for th[] */
  int              mh_count;               /* number of valid mh[] */
  QMP_mem_t       *mem[4*Q(DIM)];          /* memory for mh[] */
  int              mem_count;              /* number of valid mem[] */
  void            *send_up_buf[Q(DIM)];    /* pf up-bufs */
  void            *send_down_buf[Q(DIM)];  /* pf down-bufs */
  void            *receive_buf[2*Q(DIM)];  /* pf receive bufs (up[], down[]) */
  int              total_send;             /* bytes to send */
  int              total_receive;          /* bytes to receive */
};

/* structs ABTable and ABiTable are defined by qa0 */
struct ABTable;
struct ABiTable;

struct Q(Parameters) {
  struct Q(State) *state;
  struct ABTable  *ApTable;
  struct ABTable  *AmTable;
  struct ABTable  *AxpTable;
  struct ABTable  *AxmTable;
  struct ABTable  *BpTable;
  struct ABTable  *BmTable;
  struct ABTable  *BxpTable;
  struct ABTable  *BxmTable;
  struct ABiTable *AipTable;
  struct ABiTable *AimTable;
  struct ABiTable *BipTable;
  struct ABiTable *BimTable;
  struct ABiTable *AxipTable;
  struct ABiTable *AximTable;
  struct ABiTable *BxipTable;
  struct ABiTable *BximTable;
};

struct Q(State) {
  const char        *version;         /* to get version string into app */
  int                used;            /* gc ref counter */
  int                saved;           /* gc internal counter */
  size_t             allocated;       /* currently allocated bytes */
  size_t             max_allocated;   /* maximum allocation */

  int                error_latched;   /* if 0, allow error recording */
  int                fatal_error;     /* if 0, allow reseting latch */
  const char        *error;           /* error string */

  int                real_size;       /* 0, 4 or 8 */ 
  struct eo_lattice  even;            /* even sublattice */
  struct eo_lattice  odd;             /* odd sublattice */

  struct timeval     t0, t1;          /* for timing */
  double             time_sec;        /* seconds in the last routine */
  long long          flops;           /* FLOP in the last routine */
  long long          sent;            /* bytes sent in the last routine */
  long long          received;        /* bytes received in the last routine */

  int                Ls;              /* Ls */
  int                volume;          /* 4-d volume */
  int                lattice[Q(DIM)]; /* 4-d lattice size */
  struct local       local;           /* 4-d local sublattice */
  int                node;            /* local node id */
  int                neighbor_up[Q(DIM)];     /* the up neighbors */
  int                neighbor_down[Q(DIM)];   /* the down neighbors */
  int                network[Q(DIM)]; /* the network geometry */
  int                master_p;        /* are we the master? */
  int               *lx2v;            /* Sublattice 1-d -> 4-d translation */
  int               *v2lx;            /* Only for init */
};

typedef enum {
    CG_SUCCESS,
    CG_MAXITER,
    CG_EIGCONV,
    CG_ZEROMODE,
    CG_NOEMEM
} CG_STATUS;

/* debug printing */
extern int QDP_this_node;
extern int QDP_is_initialized(void);

/* Deflator state */
#include "deflator-la.h"

#if defined(HAVE_LAPACK)
#elif defined(HAVE_GSL)
#  include <gsl/gsl_vector.h>
#  include <gsl/gsl_matrix.h>
#  include <gsl/gsl_eigen.h>
#else
#  error "no linear algebra library"
#endif

#define DEFLATOR_VEC_SIZE(pstate) ((pstate)->even.full_size)
struct Q(Deflator) {
    struct Q(State) *state;

    int                 dim;     /* size of problem vectors, == State.size */
    int                 Ls;      /* flavor dimension extend */

    int                 vmax;
    int                 vsize;
    int                 nev;
    int                 umax;
    int                 usize;
    int                 frozen;

    /* eig current state */
    double              eps;
    double              resid_norm_sq_min;
    
    latmat_c            V;
    doublecomplex       *T;

    /* incr_eig current state */
    latmat_c            U;
    doublecomplex       *H;
    doublecomplex       *C;


    long int            lwork;
    doublecomplex       *zwork;
    doublecomplex       *hevecs2;
    double              *hevals;
#if defined(HAVE_LAPACK)
    doublecomplex       *hevecs1;
    doublecomplex       *tau;
    double              *rwork;
    
    double              *debug_hevals;
    long int            debug_lwork;
    doublecomplex       *debug_zwork;
    double              *debug_rwork;

#elif defined(HAVE_GSL)
    doublecomplex       *zwork2;
    gsl_matrix_complex  *gsl_T_full;
    gsl_matrix_complex  *gsl_hevecs1;
    gsl_vector          *gsl_hevals1;
    gsl_eigen_hermv_workspace *gsl_wkspace1;
    gsl_matrix_complex  *gsl_T_m1;
    gsl_matrix_complex  *gsl_hevecs2;
    gsl_vector          *gsl_hevals2;
    gsl_eigen_hermv_workspace *gsl_wkspace2;
    gsl_matrix_complex  *gsl_T_proj;
    gsl_matrix_complex  *gsl_hevecs3;
    gsl_eigen_hermv_workspace *gsl_wkspace3;
    gsl_matrix_complex  *gsl_QR;
    gsl_matrix_complex  *gsl_Q_unpack;
    gsl_matrix_complex  *gsl_tmp_MxS;
    gsl_vector_complex  *gsl_tau;
    size_t              *hevals_select1;
    size_t              *hevals_select2;

    gsl_vector          *debug_gsl_hevals;
    gsl_eigen_herm_workspace *debug_gsl_wkspace;

#else
#  error "no linear algebra library"
#endif

    latmat_c            tmp_V;
    latvec_c            work_c_1;
    latvec_c            work_c_2;
    latvec_c            work_c_3;

};

/* Timing */
#define BEGIN_TIMING(s) do { gettimeofday(&((s)->t0), NULL); } while (0)
#define END_TIMING(s, f, snd, rcv) do { \
    gettimeofday(&((s)->t1), NULL); \
    (s)->time_sec = ((s)->t1.tv_sec - (s)->t0.tv_sec) \
      + 1e-6 * ((s)->t1.tv_usec - (s)->t0.tv_usec); \
    (s)->flops = (f); (s)->sent = (snd); (s)->received = (rcv); } while (0)

/* Argument checking */
#define DECLARE_STATE struct Q(State) *state = NULL
#define CHECK_ARG0(n) do { if ((n) == 0) return 1;      \
    state = (n)->state; } while (0)
#define CHECK_ARGn(n,f) do { if ((n) == 0)                              \
      return q(set_error)(state, 0, f "(): NULL argument");             \
    if ((n)->state != state)                                            \
      return q(set_error)(state, 0, f "(): geometry mismatch"); } while (0)
#define CHECK_POINTER(n,f) do { if ((n) == 0)                           \
      return q(set_error)(state, 0, f "(): NULL argument"); } while (0)

#endif /* !defined(MARK_B9BA8123_0F1A_40FD_8827_42266FE32F3E) */

# undef qx
# undef QX
# undef REAL
# undef SUn
# undef Fermion
# undef VectorFermion
# undef ProjectedFermion
# undef MxM_workspace
# if QOP_MDWF_DEFAULT_PRECISION=='D'
#  define qx(x) qop_d3_mdwf_##x
#  define QX(x) QOP_D3_MDWF_##x
#  define REAL double
#  define SUn  SUnD
#  define Fermion FermionD
#  define VectorFermion VectorFermionD
#  define ProjectedFermion ProjectedFermionD
#  define MxM_workspace MxM_workspaceD
# endif
# if QOP_MDWF_DEFAULT_PRECISION=='F'
#  define qx(x) qop_f3_mdwf_##x
#  define QX(x) QOP_F3_MDWF_##x
#  define REAL float
#  define SUn  SUnF
#  define Fermion FermionF
#  define VectorFermion VectorFermionF
#  define ProjectedFermion ProjectedFermionF
#  define MxM_workspace MxM_workspaceF
# endif

/* MDWF types */
struct QX(Fermion) {
  struct Q(State) *state;
  size_t size;
  struct Fermion *even;  
  struct Fermion *odd;
};

struct QX(HalfFermion) {
  struct Q(State) *state;
  size_t size;
  struct Fermion *even;  
};

struct QX(VectorFermion) {
  struct Q(State) *state;
  size_t size;
  int    count;
  struct VectorFermion *even;  
};

struct QX(Gauge) {
  struct Q(State) *state;
  size_t size;
  struct SUn *data;
};

/* mixed precision operations */
/* Fd = Fd + Ff */
unsigned int q(f_d_peq_f)(struct FermionD *dst,
                          int size, int Ls,
                          const struct FermionF *src_f);
/* Ff = Fd - Fd */
unsigned int q(f_f_eq_dmd)(struct FermionF *dst,
                           int size, int Ls,
                           const struct FermionD *src_a,
                           const struct FermionD *src_b);

unsigned int q(f_f_eq_dmd_norm2)(struct FermionF *dst,
                                 double *local_norm2,
                                 int size, int Ls,
                                 const struct FermionD *src_a,
                                 const struct FermionD *src_b);

/* converting gauge from double down to float */
void q(g_f_eq_d)(struct SUnF *dst,
                 int size,
                 const struct SUnD *src);

/* the mixed solver */
int q(mixed_cg)(struct Q(State)             *state,
                const char                  *name,
                const struct Q(Parameters)  *params,
                struct QD(Fermion)          *psi,
                int                         *out_iterations,
                double                      *out_epsilon,
                const struct QD(Fermion)    *psi_0,
                const struct QD(Gauge)      *gauge,
                const struct QD(Fermion)    *eta,
                struct Q(Deflator)          *deflator,
                int                          f_iter,
                double                       f_epsilon,
                int                          max_iterations,
                double                       min_epsilon,
                unsigned int                 options);

/* handling eig deflator */
int Q(create_deflator)(
        struct Q(Deflator) **deflator_ptr,
        struct Q(State) *s,
        int vmax, int nev,
        double eps, int umax);
void Q(free_deflator)(struct Q(Deflator) **deflator_ptr);
void Q(deflator_reset)(struct Q(Deflator) *deflator);
void Q(deflator_stop)(struct Q(Deflator) *deflator);
void Q(deflator_resume)(struct Q(Deflator) *deflator);
int Q(deflator_is_stopped)(struct Q(Deflator) *deflator);
int q(df_preamble)(struct Q(State)           *state,
                   struct Q(Deflator)        *deflator,
                   struct FermionF           *psi_e,
                   struct FermionF           *rho_e,
                   double                    *rho_norm2,
                   struct FermionF           *chi_e, /* const ! */
                   struct MxM_workspaceF     *ws,
                   unsigned int              options);
int q(df_update0)(struct Q(State)          *state,
                  struct Q(Deflator)       *deflator,
                  double                    a1,
                  double                    b1,
                  double                    a0,
                  double                    b0,
                  double                    r,
                  struct FermionF          *rho,
                  unsigned int              options);
int q(df_update1)(struct Q(State)          *state,
                  struct Q(Deflator)       *deflator,
                  double                    a1,
                  double                    b1,
                  double                    a0,
                  double                    b0,
                  double                    r,
                  struct FermionF          *rho,
                  struct FermionF          *A_rho,
                  unsigned int              options);
int q(df_postamble)(struct Q(State)           *state,
                    struct Q(Deflator)        *deflator,
                    struct MxM_workspaceF     *ws,
                    unsigned int               options);

/* layout translation */
void q(l2v)(int x[Q(DIM)], const struct local *local, int p);
int q(v2l)(const int x[Q(DIM)], const struct local *local);

/* Implementation functions */
int q(set_error)(struct Q(State) *state, int fatal, const char *error);

int q(setup_comm)(struct Q(State) *state, int real_size);
int q(free_comm)(struct Q(State) *state);

void *q(malloc)(struct Q(State) *state, size_t bytes);
void *q(allocate_aligned)(struct Q(State) *state,
                          size_t *size, void **aligned_ptr,
                          size_t hdr_size, size_t bulk_size);
void *qx(allocate_eo)(struct Q(State) *state,
                      size_t *size, void **aligned_ptr,
                      size_t hdr_size, int even_count, int odd_count);
void q(free)(struct Q(State) *state, void *ptr, size_t bytes);
void q(cleanup_state)(struct Q(State) *state);
void *qx(step_even)(struct Q(State) *state, void *aligned_ptr);
void *qx(step_odd)(struct Q(State) *state, void *aligned_ptr);

void qx(x_import)(struct eo_lattice *eo,
                  double r[],
                  struct Fermion *data, 
                  void (*reader)(double *val_re,
                                 double *val_im,
                                 const int pos[Q(DIM)+1],
                                 int color,
                                 int dirac, 
                                 void *env),
                  void *env);
void qx(x4_import)(struct eo_lattice *eo,
                   double r[],
                   struct Fermion *data, 
                   void (*reader)(double *val_re,
                                  double *val_im,
                                  const int pos[Q(DIM)+1],
                                  int color,
                                  int dirac, 
                                  void *env),
                  void *env);
void qx(x_export)(struct eo_lattice *eo,
                  double r[],
                  const struct Fermion *data, 
                  void (*writer)(const int pos[Q(DIM)+1],
                                 int color,
                                 int dirac, 
                                 double val_re,
                                 double val_im,
                                 void *env),
                  void *env);
void qx(x4_export)(struct eo_lattice *eo,
                   double r[],
                   const struct Fermion *data, 
                   void (*writer)(const int pos[Q(DIM)],
                                  int color,
                                  int dirac, 
                                  double val_re,
                                  double val_im,
                                  void *env),
                   void *env);
void qx(x_midpoint)(struct eo_lattice *eo,
                    double r[],
                    const struct Fermion *data,
                    void (*writer)(const int pos[Q(DIM)],
                                   double value,
                                   void *env),
                    void *env);

/* Projections */
typedef unsigned int (*qx(Up_project))(struct ProjectedFermion *r,
                                       int size, int Ls,
                                       const struct up_pack *link,
                                       const struct SUn *U,
                                       const struct Fermion *f);
typedef unsigned int (*qx(Down_project))(struct ProjectedFermion *r,
                                         int size, int Ls,
                                         const struct down_pack *link,
                                         const struct Fermion *f);
unsigned int qx(proj_g0plus)(struct ProjectedFermion *r,
                             int size, int Ls,
                             const struct down_pack *link,
                             const struct Fermion *f);
unsigned int qx(proj_g1plus)(struct ProjectedFermion *r,
                             int size, int Ls,
                             const struct down_pack *link,
                             const struct Fermion *f);
unsigned int qx(proj_g2plus)(struct ProjectedFermion *r,
                             int size, int Ls,
                             const struct down_pack *link,
                             const struct Fermion *f);
unsigned int qx(proj_g3plus)(struct ProjectedFermion *r,
                             int size, int Ls,
                             const struct down_pack *link,
                             const struct Fermion *f);
unsigned int qx(proj_g0minus)(struct ProjectedFermion *r,
                              int size, int Ls,
                              const struct down_pack *link,
                              const struct Fermion *f);
unsigned int qx(proj_g1minus)(struct ProjectedFermion *r,
                              int size, int Ls,
                              const struct down_pack *link,
                              const struct Fermion *f);
unsigned int qx(proj_g2minus)(struct ProjectedFermion *r,
                              int size, int Ls,
                              const struct down_pack *link,
                              const struct Fermion *f);
unsigned int qx(proj_g3minus)(struct ProjectedFermion *r,
                              int size, int Ls,
                              const struct down_pack *link,
                              const struct Fermion *f);
unsigned int qx(proj_Ucg0plus)(struct ProjectedFermion *r,
                               int size, int Ls,
                               const struct up_pack *link,
                               const struct SUn *U,
                               const struct Fermion *f);
unsigned int qx(proj_Ucg1plus)(struct ProjectedFermion *r,
                               int size, int Ls,
                               const struct up_pack *link,
                               const struct SUn *U,
                               const struct Fermion *f);
unsigned int qx(proj_Ucg2plus)(struct ProjectedFermion *r,
                               int size, int Ls,
                               const struct up_pack *link,
                               const struct SUn *U,
                               const struct Fermion *f);
unsigned int qx(proj_Ucg3plus)(struct ProjectedFermion *r,
                               int size, int Ls,
                               const struct up_pack *link,
                               const struct SUn *U,
                               const struct Fermion *f);
unsigned int qx(proj_Ucg0minus)(struct ProjectedFermion *r,
                                int size, int Ls,
                                const struct up_pack *link,
                                const struct SUn *U,
                                const struct Fermion *f);
unsigned int qx(proj_Ucg1minus)(struct ProjectedFermion *r,
                                int size, int Ls,
                                const struct up_pack *link,
                                const struct SUn *U,
                                const struct Fermion *f);
unsigned int qx(proj_Ucg2minus)(struct ProjectedFermion *r,
                                int size, int Ls,
                                const struct up_pack *link,
                                const struct SUn *U,
                                const struct Fermion *f);
unsigned int qx(proj_Ucg3minus)(struct ProjectedFermion *r,
                                int size, int Ls,
                                const struct up_pack *link,
                                const struct SUn *U,
                                const struct Fermion *f);

/* projection tables */
/*  normal projection */
extern qx(Up_project) qx(up_project_n)[Q(DIM)];
extern qx(Down_project) qx(down_project_n)[Q(DIM)];
/*  conjugated projection */
extern qx(Up_project) qx(up_project_x)[Q(DIM)];
extern qx(Down_project) qx(down_project_x)[Q(DIM)];
/* compute projections on the boundary and fill the send buffers */
void qx(boundary)(struct eo_lattice *xy,
                  int Ls,
                  const qx(Up_project) up_proj[],
                  const qx(Down_project) down_proj[],
                  const struct SUn *U,
                  const struct Fermion *src_y,
                  long long *flops);
/* same as above by do the down boundary only */
void qx(down_boundary)(struct eo_lattice *xy,
                       int Ls,
                       const qx(Down_project) down_proj[],
                       const struct Fermion *src_y,
                       long long *flops);

/* Backend controled structure sizes */
int q(sizeof_neighbor)(int volume);
int q(sizeof_up_pack)(int volume);
int q(sizeof_down_pack)(int volume);
int q(sizeof_ABTable)(int Ls);
int q(sizeof_ABiTable)(int Ls);
int qx(sizeof_fermion)(int volume, int Ls);
int qx(sizeof_projected_fermion)(int volume, int Ls);
int qx(sizeof_gauge)(int volume);
int qx(sizeof_vfermion)(int volume, int Ls, int count);

/* qa0 level data access routines */
void qx(put_gauge)(struct SUn *ptr, int pos, const double r[]);
void qx(put_fermion)(struct Fermion *data, int pos, int Ls, const double r[]);
void qx(get_fermion)(double r[], const struct Fermion *data, int pos, int Ls);
int q(get_down_pack_f)(const struct down_pack *up, int p);
int q(get_up_pack_f)(const struct up_pack *up, int p);
void q(put_down_pack)(struct down_pack *down, int p, int f);
void q(get_down_pack)(int *f, const struct down_pack *up, int p);
void q(put_up_pack)(struct up_pack *up, int p, int f, int u);
void q(get_up_pack)(int *f, int *u, const struct up_pack *up, int p);
void q(put_neighbor)(struct neighbor *n, int p,
                     int m,
                     const int f_up[Q(DIM)], int u_up,
                     const int f_down[Q(DIM)], const int u_down[Q(DIM)]);
void q(get_neighbor)(int *m, int *f_up, int *u_up,
                     int *f_down, int *u_down,
                     const struct neighbor *n, int p);
void q(fix_neighbor_f_up)(struct neighbor *n, int p, int f_up, int d);
void q(fix_neighbor_f_down)(struct neighbor *n, int p, int f_down, int d);
void q(put_ABTable)(struct ABTable *t, int i, double w, double v);
void q(put_ABiTableZ)(struct ABiTable *t, double z);
void q(put_ABiTable)(struct ABiTable *t,
                     int i, double vp, double sp, double fp);

/* Linear algebra on fermions */
void qx(f_zero)(struct Fermion *dst, 
                int size, int Ls);
void qx(f_copy)(struct Fermion *dst, 
                int size, int Ls,
                const struct Fermion *src);
unsigned int qx(f_dot)(double *v_r, double *v_i,
                       int size, int Ls,
                       const struct Fermion *a,
                       const struct Fermion *b);
unsigned int qx(f_add3)(struct Fermion *r,
                        int size, int Ls,
                        const struct Fermion *a,
                        double s,
                        const struct Fermion *b);
unsigned int qx(f_add2)(struct Fermion *r,
                        int size, int Ls,
                        double s,
                        const struct Fermion *b);
unsigned int qx(f_cadd2)(struct Fermion *r,
                         int size, int Ls,
                         double sr, double si,
                         const struct Fermion *b);
unsigned int qx(f_add2_norm)(struct Fermion *r,
                             double *local_norm,
                             int size, int Ls,
                             double s,
                             const struct Fermion *b);
unsigned int qx(f_rmul1)(struct Fermion *r,
                         int size, int Ls,
                         double s);
unsigned int qx(f_add2x)(struct Fermion *r,
                         int size, int Ls,
                         double s,
                         const struct Fermion *b);
unsigned int qx(f_norm)(double *s,
                        int size, int Ls,
                        const struct Fermion *a);
unsigned int qx(f_diff_norm)(double *s,
                             int size, int Ls,
                             const struct Fermion *a,
                             const struct Fermion *b);
void qx(fv_zero)(struct VectorFermion *vf,
                 int size, int Ls, int count);
void qx(fv_copy)(struct VectorFermion *vf,
                 int size, int Ls, int count,
                 const struct Fermion *f);
void qx(fv_get)(struct Fermion *f,
                int size, int Ls, int count,
                const struct VectorFermion *vf, int k);
void qx(fv_put)(struct VectorFermion *vf, int k,
                int size, int Ls, int count,
                const struct Fermion *f);

/* algebra for arrays of fermions */

unsigned int qx(vf_zero)(struct vFermion *dst, 
                         int size, int Ls, int len);

/* fv[fv_begin + (0 .. len-1)] = gv[gv_begin + (0 .. len-1)]
*/
unsigned int qx(vf_copy)(int size, int Ls, int len,
                         struct vFermion *fv, int fv_size, int fv_begin,
                         const struct vFermion *gv, int gv_size, int gv_begin);
/*
 * set fv[idx] = x
*/
unsigned int qx(vf_put)(int size, int Ls,
                        struct vFermion *fv, int fv_size, int fv_idx,
                        const struct Fermion *x);

/*
 * read x = fv[idx]
*/
unsigned int qx(vf_get)(int size, int Ls,
                        struct Fermion *x,
                        const struct vFermion *fv, int fv_size, int fv_idx);

/*
*   g = fv[fv_begin + (0 .. f_vlen-1)] . v
*   v is a complex vector [fv_len] indexed as [re:0/im:1 + 2 * i]
*/
unsigned int qx(vf_dot_vz)(int size, int Ls,
                           struct Fermion *g,
                           const struct vFermion *fv,
                           int fv_size, int fv_begin, int fv_len,
                           const double *v);

/*
*   gv[gv_begin + (0 .. gv_len-1)] = fv[fv_begin + (0 .. f_len - 1)] . m
*   m is a complex matrix [fv_len*gv_len] indexed as [re:0/im:1 + 2 * (row + ldm * col) ]
*/
unsigned int qx(vf_dot_mz)(int size, int Ls,
                           struct vFermion *gv,
                           int gv_row_size, int gv_begin, int gv_len,
                           const struct vFermion *fv,
                           int fv_row_size, int fv_begin, int fv_len,
                           const double *m, int ldm);

/*  This includes global reduction
 *  c[i] = herm(fv[fv_begin+i]) * g 
 *      for all i = (0 .. fv_len-1)
 *  c is complex vector as [re:0/im:1 + 2 * i]
 */
unsigned int qx(vfH_dot_f)(int size, int Ls,
                           double *c,
                           const struct vFermion *fv,
                           int fv_size, int fv_begin, int fv_len,
                           const struct Fermion *g);

/* Local part of the above */
unsigned int qx(do_vfH_dot_f)(int size, int Ls,
                              double *c,
                              const struct vFermion *fv,
                              int fv_size, int fv_begin, int fv_len,
                              const struct Fermion *g);

/* This includes global reduction
 * c[i,j] = herm(fv[fv_begin + i]) . g[gv_begin+j] 
 *      for all i = (0 .. fv_len-1), 
 *              j = (0 .. gv_len-1),
 * c is a complex matrix as [re:0/im:1 + 2 * (i + ldc * j)]
 */
unsigned int qx(vfH_dot_vf)(int size, int Ls,
                            double *c, int ldc,
                            const struct vFermion *fv,
                            int fv_size, int fv_begin, int fv_len,
                            const struct vFermion *gv,
                            int gv_size, int gv_begin, int gv_len);

/* local part of the above */
unsigned int qx(do_vfH_dot_vf)(int size, int Ls,
                               double *c, int ldc,
                               const struct vFermion *fv,
                               int fv_size, int fv_begin, int fv_len,
                               const struct vFermion *gv,
                               int gv_size, int gv_begin, int gv_len);

/* basic matrices */
unsigned int qx(op_norm2)(double *global_norm,
                          const struct QX(Fermion) *psi,
                          struct Q(State) *state);
unsigned int qx(do_A)(struct Fermion *r_x,
                      int size, int Ls,
                      const struct ABTable *aptable,
                      const struct ABTable *amtable,
                      const struct Fermion *s_x);
unsigned int qx(do_A_conj)(struct Fermion *r_x,
                           int size, int Ls,
                           const struct ABTable *axptable,
                           const struct ABTable *axmtable,
                           const struct Fermion *s_x);
unsigned int qx(do_A_inverse)(struct Fermion *r,
                              int size, int Ls,
                              const struct ABiTable *iatable_p,
                              const struct ABiTable *iatable_m,
                              const struct Fermion *x);
unsigned int qx(do_A_conj_inverse)(struct Fermion *r,
                                   int size, int Ls,
                                   const struct ABiTable *iatable_p,
                                   const struct ABiTable *iatable_m,
                                   const struct Fermion *x);
unsigned int qx(do_F)(struct Fermion *res_x,
                      int start, int size, int Ls,
                      const struct neighbor *neighbor,
                      const struct SUn *U,
                      const struct Fermion *src_y,
                      void *rb[]);
unsigned int qx(do_F_conj)(struct Fermion *res_x,
                           int start, int size, int Ls,
                           const struct neighbor *neighbor,
                           const struct SUn *U,
                           const struct Fermion *src_y,
                           void *rb[]);

/* basic axial current */
unsigned int qx(do_axial_current)(double *val, /* [ 2 * Q(DIM) ] */
                                  int p,
                                  int Ls,
                                  const struct neighbor *neighbor,
                                  const struct SUn *U,
                                  const struct Fermion *s_x,
                                  const struct Fermion *s_y,
                                  void *rv[]);
/* basic A+F, A and B */
unsigned int qx(do_ApF)(struct Fermion *r_x,
                        int start, int size, int Ls,
                        const struct ABTable *aptable,
                        const struct ABTable *amtable,
                        const struct neighbor *neighbor,
                        const struct SUn *U,
                        const struct Fermion *s_x,
                        const struct Fermion *s_y,
                        void *rb[]);
unsigned int qx(do_ApF_norm)(struct Fermion *r_x,
                             double *local_norm,
                             int start, int size, int Ls,
                             const struct ABTable *aptable,
                             const struct ABTable *amtable,
                             const struct neighbor *neighbor,
                             const struct SUn *U,
                             const struct Fermion *s_x,
                             const struct Fermion *s_y,
                             void *rb[]);
unsigned int qx(do_AxpBxFx)(struct Fermion *r_x,
                            int start, int size, int Ls,
                            const struct ABTable *aptable,
                            const struct ABTable *amtable,
                            const struct ABTable *bptable,
                            const struct ABTable *bmtable,
                            const struct neighbor *neighbor,
                            const struct SUn *U,
                            const struct Fermion *s_x,
                            const struct Fermion *s_y,
                            void *rb[]);
/* basic B 1/A F */
unsigned int qx(do_BA1)(struct Fermion *r_x,
                        int size, int Ls,
                        const struct ABTable *bptable,
                        const struct ABTable *bmtable,
                        const struct ABiTable *iatable_p,
                        const struct ABiTable *iatable_m,
                        const struct Fermion *s_x);
unsigned int qx(do_BA1F)(struct Fermion *r_y,
                         int start, int size, int Ls,
                         const struct ABTable *bptable,
                         const struct ABTable *bmtable,
                         const struct ABiTable *iatable_p,
                         const struct ABiTable *iatable_m,
                         const struct neighbor *neighbor,
                         const struct SUn *U,
                         const struct Fermion *s_x,
                         void *rb[]);
unsigned int qx(do_1mBA1F)(struct Fermion *r_y,
                           int start, int size, int Ls,
                           const struct ABTable *bptable,
                           const struct ABTable *bmtable,
                           const struct ABiTable *iatable_p,
                           const struct ABiTable *iatable_m,
                           const struct neighbor *neighbor,
                           const struct SUn *U,
                           const struct Fermion *a_y,
                           const struct Fermion *b_x,
                           void *rb[]);
unsigned int qx(do_1mBA1F_norm)(struct Fermion *r_y,
                                double *local_norm,
                                int start, int size, int Ls,
                                const struct ABTable *bptable,
                                const struct ABTable *bmtable,
                                const struct ABiTable *iatable_p,
                                const struct ABiTable *iatable_m,
                                const struct neighbor *neighbor,
                                const struct SUn *U,
                                const struct Fermion *a_y,
                                const struct Fermion *b_x,
                                void *rb[]);
unsigned int qx(do_A1xBx)(struct Fermion *r_y,
                          int size, int Ls,
                          const struct ABTable *bptable,
                          const struct ABTable *bmtable,
                          const struct ABiTable *iatable_p,
                          const struct ABiTable *iatable_m,
                          const struct Fermion *b_y);
unsigned int qx(do_A1xBxFx)(struct Fermion *r_x,
                            int start, int size, int Ls,
                            const struct ABiTable *aiptable,
                            const struct ABiTable *aimtable,
                            const struct ABTable *bptable,
                            const struct ABTable *bmtable,
                            const struct neighbor *neighbor,
                            const struct SUn *U,
                            const struct Fermion *s_y,
                            void *rb[]);
unsigned int qx(do_1mF)(struct Fermion *r_y,
                        int start, int size, int Ls,
                        const struct neighbor *neighbor,
                        const struct SUn *U,
                        const struct Fermion *a_y,
                        const struct Fermion *b_x,
                        void *rb[]);
unsigned int qx(do_1mFx)(struct Fermion *r_y,
                         int start, int size, int Ls,
                         const struct neighbor *neighbor,
                         const struct SUn *U,
                         const struct Fermion *a_y,
                         const struct Fermion *b_x,
                         void *rb[]);
unsigned int qx(do_1mFx_norm)(struct Fermion *r_y,
                              double *local_norm,
                              int start, int size, int Ls,
                              const struct neighbor *neighbor,
                              const struct SUn *U,
                              const struct Fermion *a_y,
                              const struct Fermion *b_x,
                              void *rb[]);

/* even/odd level routines */
void qx(op_A)(struct Fermion *r_x,
              struct eo_lattice *xy,
              const struct Q(Parameters) *params,
              const struct Fermion *s_x,
              long long *flops);
void qx(op_Ax)(struct Fermion *r_x,
               struct eo_lattice *xy,
               const struct Q(Parameters) *params,
               const struct Fermion *s_x,
               long long *flops);
void qx(op_A1)(struct Fermion *r_x,
               struct eo_lattice *xy,
               const struct Q(Parameters) *params,
               const struct Fermion *s_x,
               long long *flops);
void qx(op_A1x)(struct Fermion *r_x,
                struct eo_lattice *xy,
                const struct Q(Parameters) *params,
                const struct Fermion *s_x,
                long long *flops);
void qx(op_B)(struct Fermion *r_x,
              struct eo_lattice *xy,
              const struct Q(Parameters) *params,
              const struct Fermion *s_x,
              long long *flops);
void qx(op_Bx)(struct Fermion *r_x,
               struct eo_lattice *xy,
               const struct Q(Parameters) *params,
               const struct Fermion *s_x,
               long long *flops);
void qx(op_B1)(struct Fermion *r_x,
               struct eo_lattice *xy,
               const struct Q(Parameters) *params,
               const struct Fermion *s_x,
               long long *flops);
void qx(op_B1x)(struct Fermion *r_x,
                struct eo_lattice *xy,
                const struct Q(Parameters) *params,
                const struct Fermion *s_x,
                long long *flops);
void qx(op_A1xBx)(struct Fermion *r_x,
                  struct eo_lattice *xy,
                  const struct Q(Parameters) *params,
                  const struct Fermion *s_x,
                  long long *flops);
void qx(op_BA1)(struct Fermion *r_x,
                struct eo_lattice *xy,
                const struct Q(Parameters) *params,
                const struct Fermion *a_x,
                long long *flops);
void qx(op_F)(struct Fermion *r_x,
              struct eo_lattice *xy,
              const struct SUn *U,
              const struct Fermion *s_y,
              long long *flops,
              long long *sent,
              long long *received);
void qx(op_Fx)(struct Fermion *r_x,
               struct eo_lattice *xy,
               const struct SUn *U,
               const struct Fermion *s_y,
               long long *flops,
               long long *sent,
               long long *received);
void qx(op_ApF)(struct Fermion *r_x,
                struct eo_lattice *xy,
                const struct Q(Parameters) *params,
                const struct SUn *U,
                const struct Fermion *a_x,
                const struct Fermion *a_y,
                long long *flops,
                long long *sent,
                long long *received);
void qx(op_ApF_norm)(struct Fermion *r_x,
                     double *local_norm,
                     struct eo_lattice *xy,
                     const struct Q(Parameters) *params,
                     const struct SUn *U,
                     const struct Fermion *a_x,
                     const struct Fermion *a_y,
                     long long *flops,
                     long long *sent,
                     long long *received);
void qx(op_AxpBxFx)(struct Fermion *r_x,
                    struct eo_lattice *xy,
                    const struct Q(Parameters) *params,
                    const struct SUn *U,
                    const struct Fermion *a_x,
                    const struct Fermion *a_y,
                    long long *flops,
                    long long *sent,
                    long long *received);
void qx(op_1mF)(struct Fermion *r_x,
                struct eo_lattice *xy,
                const struct SUn *U,
                const struct Fermion *a_x,
                const struct Fermion *a_y,
                long long *flops,
                long long *sent,
                long long *received);
void qx(op_1mFx)(struct Fermion *r_x,
                 struct eo_lattice *xy,
                 const struct SUn *U,
                 const struct Fermion *a_x,
                 const struct Fermion *a_y,
                 long long *flops,
                 long long *sent,
                 long long *received);
void qx(op_1mBA1F)(struct Fermion *r_x,
                   struct eo_lattice *xy,
                   const struct Q(Parameters) *params,
                   const struct SUn *U,
                   const struct Fermion *a_x,
                   const struct Fermion *a_y,
                   long long *flops,
                   long long *sent,
                   long long *received);
void qx(op_1mBA1F_norm)(struct Fermion *r_x,
                        double *norm,
                        struct eo_lattice *xy,
                        const struct Q(Parameters) *params,
                        const struct SUn *U,
                        const struct Fermion *a_x,
                        const struct Fermion *a_y,
                        long long *flops,
                        long long *sent,
                        long long *received);
void qx(op_1mFx_norm)(struct Fermion *r_x,
                      double *local_norm,
                      struct eo_lattice *xy,
                      const struct SUn *U,
                      const struct Fermion *a_x,
                      const struct Fermion *a_y,
                      long long *flops,
                      long long *sent,
                      long long *received);
void qx(op_A1xBxFx)(struct Fermion *r_x,
                    struct eo_lattice *xy,
                    const struct Q(Parameters) *params,
                    const struct SUn *U,
                    const struct Fermion *a_y,
                    long long *flops,
                    long long *sent,
                    long long *received);
void qx(op_BA1F)(struct Fermion *r_x,
                 struct eo_lattice *xy,
                 const struct Q(Parameters) *params,
                 const struct SUn *U,
                 const struct Fermion *a_y,
                 long long *flops,
                 long long *sent,
                 long long *received);
void qx(op_axial_current)(void (*writer)(const int pos[Q(DIM)],
                                         int dir,
                                         double value,
                                         void *env),
                          void *env,
                          struct eo_lattice *xy,
                          struct eo_lattice *yx,
                          const struct SUn *U,
                          const struct Fermion *a_x,
                          const struct Fermion *a_y,
                          long long *flops,
                          long long *sent,
                          long long *received,
    int node);
void qx(op_D)(struct Fermion *r_x,
              struct eo_lattice *xy,
              struct eo_lattice *yx,
              const struct Q(Parameters) *params,
              const struct SUn *U,
              const struct Fermion *a_x,
              const struct Fermion *a_y,
              long long *flops,
              long long *sent,
              long long *received,
              struct Fermion *tmp_y);
void qx(op_D_norm)(struct Fermion *r_x,
                   double *local_norm,
                   struct eo_lattice *xy,
                   struct eo_lattice *yx,
                   const struct Q(Parameters) *params,
                   const struct SUn *U,
                   const struct Fermion *a_x,
                   const struct Fermion *a_y,
                   long long *flops,
                   long long *sent,
                   long long *received,
                   struct Fermion *tmp_y);
void qx(op_even_M)(struct Fermion *r_x,
                   struct Q(State) *state,
                   const struct Q(Parameters) *params,
                   const struct SUn *U,
                   const struct Fermion *a_x,
                   long long *flops,
                   long long *sent,
                   long long *received,
                   struct Fermion *tmp_y);
void qx(op_even_Mn)(struct Fermion *r_x,
                    double *global_norm,
                    struct Q(State) *state,
                    const struct Q(Parameters) *params,
                    const struct SUn *U,
                    const struct Fermion *a_x,
                    long long *flops,
                    long long *sent,
                    long long *received,
                    struct Fermion *tmp_y);
void qx(op_even_Mx)(struct Fermion *r_x,
                    struct Q(State) *state,
                    const struct Q(Parameters) *params,
                    const struct SUn *U,
                    const struct Fermion *a_x,
                    long long *flops,
                    long long *sent,
                    long long *received,
                    struct Fermion *tmp_x,
                    struct Fermion *tmp_y);
void qx(op_even_Mxn)(struct Fermion *r_x,
                     double *global_norm,
                     struct Q(State) *state,
                     const struct Q(Parameters) *params,
                     const struct SUn *U,
                     const struct Fermion *a_x,
                     long long *flops,
                     long long *sent,
                     long long *received,
                     struct Fermion *tmp_x,
                     struct Fermion *tmp_y);
/* logging */
void qx(zprint)(struct Q(State) *state,
                const char *source,
                const char *fmt,
                ...);
/* parts of the CG solver */
void qx(cg_precondition)(struct Fermion *xi0_e,
                         struct Fermion *chi_e,
                         struct Q(State) *state,
                         const struct Q(Parameters) *params,
                         const struct SUn *U,
                         const struct Fermion *psi0_e,
                         const struct Fermion *eta_e,
                         const struct Fermion *eta_o,
                         long long *flops,
                         long long *sent,
                         long long *received,
                         struct Fermion *t0_e,
                         struct Fermion *t1_e,
                         struct Fermion *t0_o);

void qx(cg_inflate)(struct Fermion *psi_e,
                    struct Fermion *psi_o,
                    struct Q(State) *state,
                    const struct Q(Parameters) *params,
                    const struct SUn *U,
                    const struct Fermion *eta_o,
                    const struct Fermion *xi_e,
                    long long *flops,
                    long long *sent,
                    long long *received,
                    struct Fermion *t_o);
double qx(cg_dirac_error)(const struct Fermion *psi_e,
                          const struct Fermion *psi_o,
                          struct Q(State) *state,
                          const struct Q(Parameters) *params,
                          const struct SUn *U,
                          const struct Fermion *eta_e,
                          const struct Fermion *eta_o,
                          long long *flops,
                          long long *sent,
                          long long *received,
                          struct Fermion *t0_e,
                          struct Fermion *t1_e,
                          struct Fermion *t0_o);
void qx(cg_log)(double cg_res, const char *source, int iter,
                const struct Fermion *xi_e,
                struct Q(State) *state,
                const struct Q(Parameters) *params,
                const struct SUn *U,
                const struct Fermion *chi_e,
                long long *flops,
                long long *sent,
                long long *received,
                unsigned int options,
                struct Fermion *t0_e,
                struct Fermion *t1_e,
                struct Fermion *t2_e,
                struct Fermion *t0_o);

struct MxM_workspace {
    struct Q(State)             *state;
    const struct Q(Parameters)  *params;
    const struct SUn            *gauge;
    struct Fermion              *tmp_e;
    struct Fermion              *tmp2_e;
    struct Fermion              *tmp_o;
    long long                   *flops;
    long long                   *sent;
    long long                   *received;
};

void qx(cg_operator)(struct Fermion            *res_e,
                     const struct Fermion      *psi_e,
                     struct MxM_workspace      *ws);

CG_STATUS qx(cg_solver)(struct Fermion              *psi_e,
                        const char                  *source,
                        int                         *out_iter,
                        double                      *out_epsilon,
                        struct Q(State)             *state,
                        const struct Q(Parameters)  *params,
                        const struct SUn            *U,
                        const struct Fermion        *chi_e,
                        struct Q(Deflator)          *deflator,
                        int                          max_iter,
                        double                       epsilon,
                        unsigned                     options,
                        long long                   *flops,
                        long long                   *sent,
                        long long                   *received,
                        struct Fermion              *rho_e,
                        struct Fermion              *pi_e,
                        struct Fermion              *zeta_e,
                        struct Fermion              *t0_e,
                        struct Fermion              *t1_e,
                        struct Fermion              *t2_e,
                        struct Fermion              *t0_o);
int qx(scg_solver)(struct VectorFermion         *v_xi_e,
                   struct Fermion               *xi_e,
                   int                           count,
                   const char                   *source,
                   int                          *out_iterations,
                   double                       *out_epsilon,
                   struct Q(State)              *state,
                   const struct Q(Parameters)   *params,
                   const double                  shift[],
                   const struct SUn             *U,
                   const struct Fermion         *chi_e,
                   int                           max_iterations,
                   double                        min_epsilon,
                   unsigned                      options,
                   long long                    *flops,
                   long long                    *sent,
                   long long                    *received,
                   double                        dp[],
                   double                        d[],
                   double                        dn[],
                   double                        ad[],
                   double                        bdd[],
                   struct Fermion               *rho_e,
                   struct VectorFermion         *vpi_e,
                   struct Fermion               *pi_e,
                   struct Fermion               *zeta_e,
                   struct Fermion               *t0_e,
                   struct Fermion               *t1_e,
                   struct Fermion               *t2_e,
                   struct Fermion               *t0_o);
/*
 *  compute x <- x + alpha p
 *          p <- r + beta p
 */
unsigned int qx(cg_xp)(struct Fermion *x,
                       struct Fermion *p,
                       int size, int Ls,
                       double alpha,
                       double beta,
                       const struct Fermion *r);
/*
 * compute xi <- xi + a * pi
 *         v_xi[i] <- v_xi[i] + ad[i] * pi
 */
unsigned int qx(scg_madd)(struct Fermion *xi_e,
                          struct VectorFermion *v_xi_e,
                          int size, int Ls, int count,
                          double a,
                          const double *ad,
                          const struct Fermion *pi_e);
/*
 * compute xi <- xi + a * pi
 *         pi <- rho + b * pi
 *         v_xi[i] <- v_xi[i] + ad[i] * v_pi[i]
 *         v_pi[i] <- rho + bdd[i] * v_pi[i]
 */
unsigned int qx(scg_xp)(struct Fermion *xi_e,
                        struct Fermion *pi_e,
                        struct VectorFermion *v_xi_e,
                        struct VectorFermion *v_pi_e,
                        int size, int Ls, int count,
                        double a,
                        double b,
                        const double *ad,
                        const double *bdd,
                        const struct Fermion *rho_e);

/* --- other composites are here */
