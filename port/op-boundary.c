#include <mdwf.h>

void
qx(op_boundary)(struct eo_lattice *xy,
		int Ls,
		const Up_project up_proj[],
		const Down_project down_proj[],
		const struct SUn *U,
		const struct Fermion *src_y,
		long long *flops)
{
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    *flops += (up_proj[i])(xy->send_up_buf[i],
				   xy->send_up_size[i], Ls,
				   xy->up_pack[i], U, src_y);
	if (xy->send_down_size[i])
	    *flops += (down_proj[i])(xy->send_down_buf[i],
				     xy->send_down_size[i], Ls,
				     xy->down_pack[i], src_y);
    }
}
