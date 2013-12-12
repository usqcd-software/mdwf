struct QX(Fermion);
struct QX(HalfFermion);
struct QX(HalfFermionMat);
struct QX(VectorFermion);
struct QX(Gauge);
struct qx(MxM_workspace);

#include "deflator-x.h"

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
struct QX(HalfFermionMat) {
  struct Q(State) *state;
  size_t mem_size;
  qx(defl_mat) m;
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

struct qx(MxM_workspace) {
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

void *qx(allocate_eo)(struct Q(State) *state,
                      size_t *size, void **aligned_ptr,
                      size_t hdr_size, int even_count, int odd_count);
void *qx(step_even)(struct Q(State) *state, void *aligned_ptr);
void *qx(step_odd)(struct Q(State) *state, void *aligned_ptr);
void *qx(allocate_eo)(struct Q(State) *state,
                      size_t *size, void **aligned_ptr,
                      size_t hdr_size, int even_count, int odd_count);
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
int qx(sizeof_fermion)(int volume, int Ls);
int qx(sizeof_projected_fermion)(int volume, int Ls);
int qx(sizeof_gauge)(int volume);
int qx(sizeof_vfermion)(int volume, int Ls, int count);
int qx(strideof_vfermion)(int volume, int Ls);

/* qa0 level data access routines */
void qx(put_gauge)(struct SUn *ptr, int pos, const double r[]);
void qx(put_fermion)(struct Fermion *data, int pos, int Ls, const double r[]);
void qx(get_fermion)(double r[], const struct Fermion *data, int pos, int Ls);
void qx(fermion2blas)(void *data, const struct Fermion *f, int size, int Ls);
void qx(blas2fermion)(struct Fermion *f, int size, int Ls, const void *data);
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

/* fv[fv_begin + (0 .. len-1)] = gv[gv_begin + (0 .. len-1)]
*/
unsigned int qx(vf_copy)(int size, int Ls, int len,
                         struct vFermion *fv, int fv_stride, int fv_begin,
                         const struct vFermion *gv, int gv_stride, int gv_begin);
/*
 * set fv[idx] = x
*/
unsigned int qx(vf_put)(int size, int Ls,
                        struct vFermion *fv, int fv_stride, int fv_idx,
                        const struct Fermion *x);

/*
 * read x = fv[idx]
*/
unsigned int qx(vf_get)(int size, int Ls,
                        struct Fermion *x,
                        const struct vFermion *fv, int fv_stride, int fv_idx);

/*
*   g = fv[fv_begin + (0 .. f_vlen-1)] . v
*   v is a complex vector [fv_len] indexed as [re:0/im:1 + 2 * i]
*/
unsigned int qx(vf_dot_vz)(int size, int Ls,
                           struct Fermion *g,
                           const struct vFermion *fv,
                           int fv_stride, int fv_begin, int fv_len,
                           const double *v);

/*
*   gv[gv_begin + (0 .. gv_len-1)] = fv[fv_begin + (0 .. f_len - 1)] . m
*   m is a complex matrix [fv_len*gv_len] indexed as [re:0/im:1 + 2 * (row + ldm * col) ]
*/
unsigned int qx(vf_dot_mz)(int size, int Ls,
                           struct vFermion *gv,
                           int gv_stride, int gv_begin, int gv_len,
                           const struct vFermion *fv,
                           int fv_stride, int fv_begin, int fv_len,
                           const double *m, int ldm);

/*  This does not include global reduction
 *  c[i] = herm(fv[fv_begin+i]) * g 
 *      for all i = (0 .. fv_len-1)
 *  c is complex vector as [re:0/im:1 + 2 * i]
 */
unsigned int qx(do_vfH_dot_f)(int size, int Ls,
                              double *c,
                              const struct vFermion *fv,
                              int fv_stride, int fv_begin, int fv_len,
                              const struct Fermion *g);

/* This does not include global reduction
 * c[i,j] = herm(fv[fv_begin + i]) . g[gv_begin+j] 
 *      for all i = (0 .. fv_len-1), 
 *              j = (0 .. gv_len-1),
 * c is a complex matrix as [re:0/im:1 + 2 * (i + ldc * j)]
 */
unsigned int qx(do_vfH_dot_vf)(int size, int Ls,
                               double *c, int ldc,
                               const struct vFermion *fv,
                               int fv_stride, int fv_begin, int fv_len,
                               const struct vFermion *gv,
                               int gv_stride, int gv_begin, int gv_len);

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

void qx(cg_operator)(struct Fermion            *res_e,
                     const struct Fermion      *psi_e,
                     struct qx(MxM_workspace)  *ws);

CG_STATUS qx(cg_solver)(struct Fermion              *psi_e,
                        const char                  *source,
                        int                         *out_iter,
                        double                      *out_epsilon,
                        struct Q(State)             *state,
                        const struct Q(Parameters)  *params,
                        const struct SUn            *U,
                        const struct Fermion        *chi_e,
                        struct QX(Deflator)         *deflator,
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
                   double                        v[],
                   double                        w[],
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
