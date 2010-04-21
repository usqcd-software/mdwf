#include <mdwf.h>


void
qx(x4_import)(struct eo_lattice *eo,
              double r[],
              struct Fermion *data,
              void (*reader)(double *val_re,
                             double *val_im,
                             const int pos[Q(DIM)],
                             int color,
                             int dirac,
                             void *env),
              void *env)
{
  int size = eo->full_size;
  int Ls = eo->Ls;
  int p, c, d;
  int x[Q(DIM)];
#define Xidx(x4,c,d)  &(r[2*((d)+ Q(FERMION_DIM)*((c)+Q(COLORS)*(x4)))])
  double *v;

  for (p = 0; p < size; p++) {
    q(l2v)(x, eo->local, eo->lx2v[p]);
    for (c = 0; c < Q(COLORS); c++) {
        for (d = 0; d < Q(FERMION_DIM)/2; d++) {
            v = Xidx(0, c, d);
            reader(&v[0], &v[1], x, c, d, env);
        }
        for (d = Q(FERMION_DIM)/2; d < Q(FERMION_DIM); d++) {
            v = Xidx(Ls - 1, c, d);
            reader(&v[0], &v[1], x, c, d, env);
        }
    }
    qx(put_fermion)(data, p, Ls, r);
  }
}
