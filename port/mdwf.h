#ifndef MARK_B9BA8123_0F1A_40FD_8827_42266FE32F3E
#define MARK_B9BA8123_0F1A_40FD_8827_42266FE32F3E

# ifndef QOP_MDWF_DEFAULT_PRECISION
#  define QOP_MDWF_DEFAULT_PRECISION 'D'
# endif

# include <qop-mdwf3.h>
# include <stdlib.h>
# include <string.h>
# include <qmp.h>
# include <sys/time.h>

# define q(x) qop_mdwf_##x
# define Q(x) QOP_MDWF_##x
# if QOP_MDWF_DEFAULT_PRECISION=='D'
#  define qx(x) qop_d3_mdwf_##x
#  define QX(x) QOP_D3_MDWF_##x
#  define REAL double
# endif
# if QOP_MDWF_DEFAULT_PRECISION=='F'
#  define qx(x) qop_f3_mdwf_##x
#  define QX(x) QOP_F3_MDWF_##x
#  define REAL float
# endif

/* QCD types (qa0 controls these definitions) */
struct SUn;
struct Fermion;
struct ProjectedFermion;

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
  int              face_size;              /* 4-d size of the face */
  int              body_size;              /* 4-d size of the body */
  int              full_size;              /* face + body */
  int             *lx2v;                   /* 4-d layout 2 vector translation */
  int             *v2lx;                   /* only for init */
  int              Ls;                     /* Ls */
  struct local    *local;                  /* points to state.local */

  struct neighbor *body_neighbor;          /* neighbor data for the body */
  struct neighbor *face_neighbor;          /* neighbor data for the face */
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

struct QX(Gauge) {
  struct Q(State) *state;
  size_t size;
  struct SUn *data;
};

/* structs ABTable and ABiTable are defined by qa0 */
struct ABTable;
struct ABiTable;

struct Q(Parameters) {
  struct Q(State) *state;
  struct ABTable  *ATable;
  struct ABTable  *AxTable;
  struct ABTable  *BTable;
  struct ABTable  *BxTable;
  struct ABiTable *AipTable; /* XXX */
  struct ABiTable *AimTable; /* XXX */
  struct ABiTable *BipTable; /* XXX */
  struct ABiTable *BimTable; /* XXX */
};

struct Q(State) {
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
  int                node[Q(DIM)];    /* local node address */
  int                network[Q(DIM)]; /* the network geometry */
  int                master_p;        /* are we the master? */
  int               *lx2v;            /* Sublattice 1-d -> 4-d translation */
  int               *v2lx;            /* Only for init */
};

/* layout translation */
void q(l2v)(int x[Q(DIM)], const struct local *local, int p);
int q(v2l)(const int x[Q(DIM)], const struct local *local);

/* Implementation functions */
int q(setup_comm)(struct Q(State) *state, int real_size);
int q(free_comm)(struct Q(State) *state);

int q(set_error)(struct Q(State) *state, int fatal, const char *error);

void *q(malloc)(struct Q(State) *state, size_t bytes);
void *q(allocate_aligned)(struct Q(State) *state,
 			  size_t *size, void **aligned_ptr,
			  size_t hdr_size, size_t bulk_size);
void *q(allocate_eo)(struct Q(State) *state,
 	             size_t *size, void **aligned_ptr,
		     size_t hdr_size, int even_count, int odd_count,
		     size_t f_size);
void q(free)(struct Q(State) *state, void *ptr, size_t bytes);
void q(cleanup_state)(struct Q(State) *state);
void *q(step_even)(struct Q(State) *state, void *aligned_ptr, size_t fsize);
void *q(step_odd)(struct Q(State) *state, void *aligned_ptr, size_t fsize);

void q(x_import)(struct eo_lattice *eo,
		 double r[],
		 struct Fermion *data, 
		 double (*reader)(const int pos[Q(DIM)+1],
				  int color,
				  int dirac, 
				  int re_im,
				  void *env),
		 void *env);

void q(x_export)(struct eo_lattice *eo,
		 double r[],
		 const struct Fermion *data, 
		 void (*writer)(const int pos[Q(DIM)+1],
				int color,
				int dirac, 
				int re_im,
				double value,
				void *env),
		 void *env);

/* Projections */
typedef unsigned int (*Up_project)(struct ProjectedFermion *r,
				   int size, int Ls,
				   const struct up_pack *link,
				   const struct SUn *U,
				   const struct Fermion *f);
typedef unsigned int (*Down_project)(struct ProjectedFermion *r,
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

/* A+F, A and B */
unsigned int qx(do_ApF)(struct Fermion *r_x,
			int start, int size, int Ls,
			const struct ABTable *atable,
			const struct neighbor *neighbor,
			const struct SUn *U,
			const struct Fermion *s_x,
			const struct Fermion *s_y,
			void *rb[]);
#if 0
unsigned int qx(do_AxpFx)(struct Fermion *r_x,
			  int start, int size, int Ls,
			  const struct ABTable *atable,
			  const struct neighbor *neighbor,
			  const struct SUn *U,
			  const struct Fermion *s_x,
			  const struct Fermion *s_y,
			  void *rb[]);
#endif
unsigned int qx(do_A)(struct Fermion *r_x,
		      int size, int Ls,
		      const struct ABTable *atable,
		      const struct Fermion *s_x);
#if 0
unsigned int qx(do_AxpBxFx)(struct Fermion *r_x,
			    int start, int size, int Ls,
			    const struct ABTable *atable,
			    const struct ABTable *btable,
			    const struct neighbor *neighbor,
			    const struct SUn *U,
			    const struct Fermion *s_x,
			    const struct Fermion *s_y,
                            struct Fermion *tmp,
			    void *rb[]);
#endif
/* XXX  other functions */

/* Back end functions */
void qx(put_gauge)(struct SUn *ptr, int pos, const double r[]);
void qx(put_fermion)(struct Fermion *data, int pos, int Ls, const double r[]);
void qx(get_fermion)(double r[], const struct Fermion *data, int pos, int Ls);
unsigned int qx(madd_fermion)(struct Fermion *r,
			      int size, int Ls,
			      const struct Fermion *a,
			      double s,
			      const struct Fermion *b);
unsigned int qx(dot_fermion)(double *v_r, double *v_i,
			     int size, int Ls,
			     const struct Fermion *a,
			     const struct Fermion *b);
unsigned int qx(norm2_fermion)(double *v_r,
			       int size, int Ls,
			       const struct Fermion *a);
unsigned int qx(do_A_inverse)(struct Fermion *r,
			      int size, int Ls,
			      const struct ABiTable *iatable_p,
			      const struct ABiTable *iatable_m,
			      const struct Fermion *x);
/* XXX functions for cg, need better integration with the rest */
unsigned int qx(do_1AcBc)(struct Fermion *r,
			  int size, int Ls,
			  const struct ABiTable *iatable_p,
			  const struct ABiTable *iatable_m,
			  const struct ABTable *btable,
			  const struct Fermion *x);
unsigned int qx(do_B1A)(struct Fermion *r,
			int size, int Ls,
			const struct ABTable *btable,
			const struct ABiTable *iatable_p,
			const struct ABiTable *iatable_m,
			const struct Fermion *x);
unsigned int qx(do_1mF)(struct Fermion *r_x,
			int start, int size, int Ls,
			const struct neighbor *neighbor,
			const struct Fermion *s_x,
			const struct SUn *U,
			const struct Fermion *s_y,
			void *rb[]);
unsigned int qx(do_1mFc)(struct Fermion *r_y,
			 int start, int size, int Ls,
			 const struct neighbor *neighbor,
			 const struct Fermion *s_y,
			 const struct SUn *U,
			 const struct Fermion *s_x,
			 void *rb[]);
unsigned int qx(do_1AcBcFc)(struct Fermion *r_y,
			    int start, int size, int Ls,
			    const struct ABiTable *iatable_p,
			    const struct ABiTable *iatable_m,
			    const struct ABTable *btable,
			    const struct neighbor *neighbor,
			    const struct SUn *U,
			    const struct Fermion *s_x,
			    void *rb[]);
unsigned int qx(do_B1AF)(struct Fermion *r_y,
			 int start, int size, int Ls,
			 const struct ABTable *btable,
			 const struct ABiTable *iatable_p,
			 const struct ABiTable *iatable_m,
			 const struct neighbor *neighbor,
			 const struct SUn *U,
			 const struct Fermion *s_x,
			 void *rb[]);
unsigned int qx(do_1mB1AF)(struct Fermion *r_y,
			   int start, int size, int Ls,
			   const struct ABTable *btable,
			   const struct ABiTable *iatable_p,
			   const struct ABiTable *iatable_m,
			   const struct neighbor *neighbor,
			   const struct Fermion *a_y,
			   const struct SUn *U,
			   const struct Fermion *b_x,
			   void *rb[]);
unsigned int qx(do_1mB1AF_norm)(struct Fermion *r_y,
				double *norm,
				int start, int size, int Ls,
				const struct ABTable *btable,
				const struct ABiTable *iatable_p,
				const struct ABiTable *iatable_m,
				const struct neighbor *neighbor,
				const struct Fermion *a_y,
				const struct SUn *U,
				const struct Fermion *b_x,
				void *rb[]);
void qx(f_copy)(int size, int Ls,
		struct Fermion *dst, 
		const struct Fermion *src);
unsigned int qx(f_add3)(int size, int Ls,
			struct Fermion *r,
			const struct Fermion *a,
			double s,
			const struct Fermion *b);
unsigned int qx(f_add2)(int size, int Ls,
			struct Fermion *r,
			double s,
			const struct Fermion *b);
unsigned int qx(f_add2_norm)(int size, int Ls,
			     struct Fermion *r,
			     double *norm,
			     double s,
			     const struct Fermion *b);
unsigned int qx(f_add2x)(int size, int Ls,
			struct Fermion *r,
			double s,
			const struct Fermion *b);
unsigned int qx(f_norm)(int size, int Ls,
			double *s,
			const struct Fermion *a);
/* end of functions for cg */
/* Timing */
#define BEGIN_TIMING(s) do { gettimeofday(&((s)->t0), NULL); } while (0)
#define END_TIMING(s, f, snd, rcv) do { \
    gettimeofday(&((s)->t1), NULL); \
    (s)->time_sec = ((s)->t1.tv_sec - (s)->t0.tv_sec) \
      + 1e-6 * ((s)->t1.tv_usec - (s)->t0.tv_usec); \
    (s)->flops = (f); (s)->sent = (snd); (s)->received = (rcv); } while (0)

/* Argument checking */
#define DECLARE_STATE struct Q(State) *state = NULL
#define CHECK_ARG0(n) do { if ((n) == 0) return 1;	\
    state = (n)->state; } while (0)
#define CHECK_ARGn(n,f) do { if ((n) == 0)				\
      return q(set_error)(state, 0, f "(): NULL argument");		\
    if ((n)->state != state)						\
      return q(set_error)(state, 0, f "(): geometry mismatch"); } while (0)
#define CHECK_POINTER(n,f) do { if ((n) == 0)				\
      return q(set_error)(state, 0, f "(): NULL argument"); } while (0)

/* Cache size */
#define CACHE_LINE_SIZE 128
#define ALIGN(p) ((void *)((((ptrdiff_t)(p))+CACHE_LINE_SIZE-1) & \
                           ~(CACHE_LINE_SIZE-1)))

/* Backend controled structure sizes */
int q(sizeof_neighbor)(int volume);
int q(sizeof_up_pack)(int volume);
int q(sizeof_down_pack)(int volume);
int q(sizeof_ABTable)(int Ls);
int q(sizeof_ABiTable)(int Ls);
int qx(sizeof_fermion)(int volume, int Ls);
int qx(sizeof_projected_fermion)(int volume, int Ls);
int qx(sizeof_gauge)(int volume);

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
void q(put_ABTable)(struct ABTable *t,
		    int i, int ip, int im,
		    double v, double vp, double vm);
void q(put_ABiTableZ)(struct ABiTable *t, double z);
void q(put_ABiTable)(struct ABiTable *t,
		     int i,
		     double vp, double sp, double fp);

#endif /* !defined(MARK_B9BA8123_0F1A_40FD_8827_42266FE32F3E) */
