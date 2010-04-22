#include <mdwf.h>

void
qx(x_midpoint)(struct eo_lattice *eo,
               double r[],
               const struct Fermion *data,
               void (*writer)(const int pos[Q(DIM)],
                              double value,
                              void *env),
               void *env)
{
  int size = eo->full_size;
  int Ls = eo->Ls;
  int p, c, d;
  int x[Q(DIM)];
#define Xidx(x4,c,d)  &(r[2*((d)+ Q(FERMION_DIM)*((c)+Q(COLORS)*(x4)))])
  double *v;
  double s;
  int Ls2 = Ls / 2;

  for (p = 0; p < size; p++) {
    q(l2v)(x, eo->local, eo->lx2v[p]);
    qx(get_fermion)(r, data, p, Ls);
    s = 0;
    for (c = 0; c < Q(COLORS); c++) {
        for (d = 0; d < Q(FERMION_DIM) / 2; d++) {
            v = Xidx(Ls2-1, c, d);
            s += v[0] * v[0] + v[1] * v[1];
        }
        for (d = Q(FERMION_DIM) / 2; d < Q(FERMION_DIM); d++) {
            v = Xidx(Ls2, c, d);
            s += v[0] * v[0] + v[1] * v[1];
        }
    }
    writer(x, s, env);
  }
}
