#include <mdwf.h>

static unsigned int
do_ac(void                   (*writer)(const int pos[Q(DIM)],
                                       int dir,
                                       double val_re,
                                       double val_im,
                                       void *env),
      void                    *env,
      int                      start,
      int                      size,
      int                      Ls,
      struct eo_lattice       *xy,
      const struct SUn        *U,
      const struct Fermion    *s_x,
      const struct Fermion    *s_y,
      void                    *rb[],
      int                      node)
{
    int p, mu;
    double val[2 * Q(DIM)];
    int flops = 0;

    for (p = 0; p < size; p++) {
        flops += qx(do_axial_current)(val, p + start, Ls,
                                      xy->neighbor, U, s_x, s_y, rb);
        int x[Q(DIM)];
        q(l2v)(x, xy->local, xy->lx2v[p + start]);
        for (mu = 0; mu < Q(DIM); mu++) {
            writer(x, mu, val[2 * mu], val[2 * mu + 1], env);
        }
    }
    return flops;
}
                     

void
qx(op_axial_current)(void                   (*writer)(const int pos[Q(DIM)],
                                                      int dir,
                                                      double val_re,
                                                      double val_im,
                                                      void *env),
                     void                    *env,
                     struct eo_lattice       *xy,
                     struct eo_lattice       *yx,
                     const struct SUn        *U,
                     const struct Fermion    *a_x,
                     const struct Fermion    *a_y,
                     long long               *flops,
                     long long               *sent,
                     long long               *received,
                     int                      node)
{
    int Ls = xy->Ls;
    
    qx(down_boundary)(xy, Ls, qx(down_project_n), a_y, flops);

    if (xy->h_valid)
        QMP_start(xy->handle);

    *flops += do_ac(writer, env, 0, xy->body_size, Ls,
                    xy, U, a_x, a_y, NULL, node);

    if (xy->h_valid)
        QMP_wait(xy->handle);

    *flops += do_ac(writer, env, xy->body_size, xy->face_size, Ls,
                    xy, U, a_x, a_y, xy->receive_buf, node);
    *sent += xy->total_send;
    *received += xy->total_receive;
}


