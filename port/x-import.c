#include <mdwf.h>

void
qx(x_import)(struct eo_lattice *eo,
	     double r[],
	     struct Fermion *data,
	     double (*reader)(const int pos[Q(DIM)+1],
			      int color,
			      int dirac,
			      int re_im,
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
    for (v = r, x[Q(DIM)] = 0; x[Q(DIM)] < Ls; x[Q(DIM)]++) {
      for (c = 0; c < Q(COLORS); c++) {
	for (d = 0; d < Q(FERMION_DIM); d++) {
	  *v++ = reader(x, c, d, 0, env);
	  *v++ = reader(x, c, d, 1, env);
	}
      }
    }
    qx(put_fermion)(data, p, Ls, r);
  }
}
