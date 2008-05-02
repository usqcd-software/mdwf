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

/* Cache size */
#define CACHE_LINE_SIZE 128
#define ALIGN(p) ((void *)((((ptrdiff_t)(p))+CACHE_LINE_SIZE-1) & \
                           ~(CACHE_LINE_SIZE-1)))

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
int q(set_error)(struct Q(State) *state, int fatal, const char *error);

int q(setup_comm)(struct Q(State) *state, int real_size);
int q(free_comm)(struct Q(State) *state);

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

/* projection tables */
/*  normal projection */
extern Up_project qx(up_project_n)[Q(DIM)];
extern Down_project qx(down_project_n)[Q(DIM)];
/*  conjugated projection */
extern Up_project qx(up_project_x)[Q(DIM)];
extern Down_project qx(down_project_x)[Q(DIM)];
/* compute projections on the boundary and fill the send buffers */
void qx(boundary)(struct eo_lattice *xy,
		  int Ls,
		  const Up_project up_proj[],
		  const Down_project down_proj[],
		  const struct SUn *U,
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
unsigned int qx(f_add2_norm)(struct Fermion *r,
			     double *local_norm,
			     int size, int Ls,
			     double s,
			     const struct Fermion *b);
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
		const struct Fermion *eta_e,
		const struct Fermion *eta_o,
		long long *flops,
		long long *sent,
		long long *received,
		unsigned int options,
		struct Fermion *t0_e,
		struct Fermion *t1_e,
		struct Fermion *t2_e,
		struct Fermion *t0_o,
		struct Fermion *t1_o);
int qx(cg_solver)(struct Fermion *xi_e,
		  const char *source,
		  int *out_iter,
		  double *out_epsilon,
		  struct Q(State) *state,
		  const struct Q(Parameters) *params,
		  const struct SUn *U,
		  const struct Fermion *chi_e,
		  const struct Fermion *eta_e,
		  const struct Fermion *eta_o,
		  int max_iter,
		  double epsilon,
		  unsigned options,
		  long long *flops,
		  long long *sent,
		  long long *received,
		  struct Fermion *rho_e,
		  struct Fermion *pi_e,
		  struct Fermion *zeta_e,
		  struct Fermion *t0_e,
		  struct Fermion *t1_e,
		  struct Fermion *t2_e,
		  struct Fermion *t0_o,
		  struct Fermion *t1_o);
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


/* --- other composites are here */

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

#endif /* !defined(MARK_B9BA8123_0F1A_40FD_8827_42266FE32F3E) */
