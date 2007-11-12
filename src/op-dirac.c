#include <mdwf.h>

/* Conservative version with contained communications */
static Up_project up_project[Q(DIM)] = {
  qx(proj_Ucg0plus),
  qx(proj_Ucg1plus),
  qx(proj_Ucg2plus),
  qx(proj_Ucg3plus)
};

static Down_project down_project[Q(DIM)] = {
  qx(proj_g0minus),
  qx(proj_g1minus),
  qx(proj_g2minus),
  qx(proj_g3minus)
};

static void
qx(compute_ApFB)(struct Q(State) *state,
		 struct eo_lattice *xy,
		 struct eo_lattice *yx,
		 const struct Q(Parameters) *params,
		 struct Fermion *r_x,
		 const struct SUn *U,
		 const struct Fermion *s_x,
		 const struct Fermion *s_y,
		 struct Fermion *By,
		 long long *flops,
		 long long *sent,
		 long long *received)
{
  int Ls = state->Ls;
  int i;

  *flops += qx(do_A)(By, yx->full_size, Ls, params->BTable, s_y);

  for (i = 0; i < Q(DIM); i++) {
    if (xy->send_up_size[i])
      *flops += (up_project[i])(xy->send_up_buf[i],
				xy->send_up_size[i], Ls,
				xy->up_pack[i], U, By);
    if (xy->send_down_size[i])
      *flops += (down_project[i])(xy->send_down_buf[i],
				  xy->send_down_size[i], Ls,
				  xy->down_pack[i], By);
  }

  if (xy->h_valid)
    QMP_start(xy->handle);
  
  *flops += qx(do_ApF)(r_x, 0, xy->body_size, Ls,
		       params->ATable, xy->body_neighbor,
		       U, s_x, By, NULL);

  if (xy->h_valid)
    QMP_wait(xy->handle);

  *flops += qx(do_ApF)(r_x, xy->body_size, xy->face_size, Ls,
		       params->ATable, xy->face_neighbor,
		       U, s_x, By, xy->receive_buf);

  *sent += xy->total_send;
  *received += xy->total_receive;
}

int
QX(DDW_operator)(struct QX(Fermion) *result,
		 const struct Q(Parameters) *params,
		 const struct QX(Gauge) *gauge,
		 const struct QX(Fermion) *fermion)
{
  DECLARE_STATE;
  long long flops = 0;
  long long sent = 0;
  long long received = 0;
  size_t s = 0;
  void *t = NULL;
  void *ptr = NULL;
  int e_size;
  int o_size;
  int size;
  int f_size;

  CHECK_ARG0(result);
  CHECK_ARGn(params, "DDW_operator");
  CHECK_ARGn(gauge, "DDW_operator");
  CHECK_ARGn(fermion, "DDW_operator");

  if (q(setup_comm)(state, sizeof (REAL)))
    return q(set_error)(state, 0, "DDW_operator(): communication setup failed");

  e_size = state->even.full_size;
  o_size = state->odd.full_size;
  size = e_size > o_size? e_size: o_size;
  f_size = qx(sizeof_fermion)(size, state->Ls);

  ptr = q(allocate_aligned)(state, &s, &t, 0, f_size);
  if (ptr == 0)
    return q(set_error)(state, 0, "DDW_operator(): not enough memory");

  BEGIN_TIMING(state);
  qx(compute_ApFB)(state, &state->even, &state->odd, params,
		   result->even, gauge->data, fermion->even, fermion->odd, t,
		   &flops, &sent, &received);
  qx(compute_ApFB)(state, &state->odd, &state->odd, params,
		   result->odd, gauge->data, fermion->odd, fermion->even, t,
		   &flops, &sent, &received);
  END_TIMING(state, flops, sent, received);

  q(free)(state, ptr, s);

  return 0;
}	
