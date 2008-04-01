#include <mdwf.h>

/* Conservative version with contained communications */
static Up_project up_project[Q(DIM)] = {
  qx(proj_Ucg0minus),
  qx(proj_Ucg1minus),
  qx(proj_Ucg2minus),
  qx(proj_Ucg3minus)
};

static Down_project down_project[Q(DIM)] = {
  qx(proj_g0plus),
  qx(proj_g1plus),
  qx(proj_g2plus),
  qx(proj_g3plus)
};

static void
qx(compute_AxpBxFx)(struct Q(State) *state,
		    struct eo_lattice *xy,
		    struct eo_lattice *yx,
		    const struct Q(Parameters) *params,
		    struct Fermion *r_x,
		    const struct SUn *U,
		    const struct Fermion *s_x,
		    const struct Fermion *s_y,
		    struct Fermion *tmp,
		    long long *flops,
		    long long *sent,
		    long long *received)
{
  int Ls = state->Ls;
  int i;

  for (i = 0; i < Q(DIM); i++) {
    if (xy->send_up_size[i])
      *flops += (up_project[i])(xy->send_up_buf[i],
				xy->send_up_size[i], Ls,
				xy->up_pack[i], U, s_y);
    if (xy->send_down_size[i])
      *flops += (down_project[i])(xy->send_down_buf[i],
				  xy->send_down_size[i], Ls,
				  xy->down_pack[i], s_y);
  }

  if (xy->h_valid)
    QMP_start(xy->handle);
  
  *flops += qx(do_AxpBxFx)(r_x, 0, xy->body_size, Ls,
			   params->AxTable,
			   params->BxTable,
			   xy->body_neighbor,
			   U, s_x, s_y, tmp, NULL);

  if (xy->h_valid)
    QMP_wait(xy->handle);

  *flops += qx(do_AxpBxFx)(r_x, xy->body_size, xy->face_size, Ls,
			   params->AxTable,
			   params->BxTable,
			   xy->face_neighbor,
			   U, s_x, s_y, tmp, xy->receive_buf);

  *sent += xy->total_send;
  *received += xy->total_receive;
}

int
QX(DDW_operator_conjugated)(struct QX(Fermion) *result,
			    const struct Q(Parameters) *params,
			    const struct QX(Gauge) *gauge,
			    const struct QX(Fermion) *fermion)
{
    DECLARE_STATE;
    int f_size;
    size_t s;
    void *tmp = NULL;
    void *ptr = NULL;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    
    CHECK_ARG0(result);
    CHECK_ARGn(params, "DDW_operator");
    CHECK_ARGn(gauge, "DDW_operator");
    CHECK_ARGn(fermion, "DDW_operator");
    
    if (q(setup_comm)(state, sizeof (REAL)))
	return q(set_error)(state, 0, "DDW_operator_conjugated(): communication setup failed");
    
    f_size = qx(sizeof_fermion)(1, state->Ls);
    ptr = q(allocate_aligned)(state, &s, &tmp, 0, f_size);
    if (ptr == 0)
	return q(set_error)(state, 0, "DDW_operator_conjugated(): not enough memory");

    BEGIN_TIMING(state);
    qx(compute_AxpBxFx)(state, &state->even, &state->odd, params,
			result->even, gauge->data,
			fermion->even, fermion->odd, tmp,
			&flops, &sent, &received);
    qx(compute_AxpBxFx)(state, &state->odd, &state->odd, params,
			result->odd, gauge->data,
			fermion->odd, fermion->even, tmp,
			&flops, &sent, &received);
    END_TIMING(state, flops, sent, received);
    
    q(free)(state, ptr, s);

    return 0;
}	
