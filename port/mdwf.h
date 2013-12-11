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

/* ugly but working solution */
#if   QOP_MDWF_DEFAULT_PRECISION=='F'
# define QOP_MDWF_DEFAULT_PRECISION_SAVE__ 'F'
#elif QOP_MDWF_DEFAULT_PRECISION=='D'
# define QOP_MDWF_DEFAULT_PRECISION_SAVE__ 'D'
#endif

# undef QOP_MDWF_DEFAULT_PRECISION
# undef qx
# undef QX
# undef REAL
# undef SUn
# undef Fermion
# undef VectorFermion
# undef ProjectedFermion
# define QOP_MDWF_DEFAULT_PRECISION 'D'
# define qx(x) qop_d3_mdwf_##x
# define QX(x) QOP_D3_MDWF_##x
# define REAL double
# define SUn  SUnD
# define Fermion FermionD
# define VectorFermion VectorFermionD
# define ProjectedFermion ProjectedFermionD
# include "mdwf-x.h"


# undef QOP_MDWF_DEFAULT_PRECISION
# undef qx
# undef QX
# undef REAL
# undef SUn
# undef Fermion
# undef VectorFermion
# undef ProjectedFermion
# define QOP_MDWF_DEFAULT_PRECISION 'F'
# define qx(x) qop_f3_mdwf_##x
# define QX(x) QOP_F3_MDWF_##x
# define REAL float
# define SUn  SUnF
# define Fermion FermionF
# define VectorFermion VectorFermionF
# define ProjectedFermion ProjectedFermionF
# include "mdwf-x.h"


# undef QOP_MDWF_DEFAULT_PRECISION
# undef qx
# undef QX
# undef REAL
# undef SUn
# undef Fermion
# undef VectorFermion
# undef ProjectedFermion
#if   QOP_MDWF_DEFAULT_PRECISION_SAVE__=='D'
# define QOP_MDWF_DEFAULT_PRECISION 'D'
# define qx(x) qop_d3_mdwf_##x
# define QX(x) QOP_D3_MDWF_##x
# define REAL double
# define SUn  SUnD
# define Fermion FermionD
# define VectorFermion VectorFermionD
# define ProjectedFermion ProjectedFermionD
#elif QOP_MDWF_DEFAULT_PRECISION_SAVE__=='F'
# define QOP_MDWF_DEFAULT_PRECISION 'F'
# define qx(x) qop_f3_mdwf_##x
# define QX(x) QOP_F3_MDWF_##x
# define REAL float
# define SUn  SUnF
# define Fermion FermionF
# define VectorFermion VectorFermionF
# define ProjectedFermion ProjectedFermionF
#endif




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
                struct QF(Deflator)         *deflator,
                int                          f_iter,
                double                       f_epsilon,
                int                          max_iterations,
                double                       min_epsilon,
                unsigned int                 options);


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
void q(free)(struct Q(State) *state, void *ptr, size_t bytes);
void q(cleanup_state)(struct Q(State) *state);

/* Backend controled structure sizes */
int q(sizeof_neighbor)(int volume);
int q(sizeof_up_pack)(int volume);
int q(sizeof_down_pack)(int volume);
int q(sizeof_ABTable)(int Ls);
int q(sizeof_ABiTable)(int Ls);
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


/* --- other composites are here */
