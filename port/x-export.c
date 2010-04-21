#include <mdwf.h>

void
qx(x_export)(struct eo_lattice *eo,
             double r[],
             const struct Fermion *data,
             void (*writer)(const int pos[Q(DIM)+1],
                            int color,
                            int dirac,
                            double val_re,
                            double val_im,
                            void *env),
             void *env)
{
  int size = eo->full_size;
  int Ls = eo->Ls;
  int p, c, d;
  int x[Q(DIM)+1];
  double *v;

  for (p = 0; p < size; p++) {
    q(l2v)(x, eo->local, eo->lx2v[p]);
    qx(get_fermion)(r, data, p, Ls);
    for (v = r, x[Q(DIM)] = 0; x[Q(DIM)] < Ls; x[Q(DIM)]++) {
      for (c = 0; c < Q(COLORS); c++) {
          for (d = 0; d < Q(FERMION_DIM); d++, v += 2) {
              writer(x, c, d, v[0], v[1], env);
        }
      }
    }
  }
}
