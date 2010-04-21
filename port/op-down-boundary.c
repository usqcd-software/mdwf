#include <mdwf.h>

void
qx(down_boundary)(struct eo_lattice *xy,
                  int Ls,
                  const qx(Down_project) down_proj[],
                  const struct Fermion *src_y,
                  long long *flops)
{
    int i;

    for (i = 0; i < Q(DIM); i++) {
        if (xy->send_down_size[i])
            *flops += (down_proj[i])(xy->send_down_buf[i],
                                     xy->send_down_size[i], Ls,
                                     xy->down_pack[i], src_y);
    }
}
