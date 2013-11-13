#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

int
QX(deflator_add_vector)(const struct Q(Parameters)  *params,
                        const struct QX(Gauge)      *gauge,
                        struct Q(Deflator) *deflator,
                        const struct QX(HalfFermion) *hfermion)
{
  DECLARE_STATE;
  long long flops = 0;
  long long sent = 0;
  long long received = 0;
  struct MxM_workspace  ws;
  void *ptr = 0;
  void *temps = 0;
  size_t ptr_size = 0;
  latvec_c lv;

  /* check arguments */
  CHECK_ARG0(params);
  CHECK_ARGn(gauge, "deflator_add_vector");
  CHECK_ARGn(deflator, "deflator_add_vector");
  CHECK_ARGn(hfermion, "deflator_add_vector");

  if (deflator->loading == 0)
    return q(set_error)(state, 0, "deflator_add_vector(): not in loading state");

  /* allocate temps */
  ptr = qx(allocate_eo)(state, &ptr_size, &temps,
                        0,  /* header */
                        2,  /* evens */
                        1); /* odds */
  if (ptr == 0)
    return q(set_error)(state, 0, "deflator_add_vector(): not enough memory");

  BEGIN_TIMING(state);

  ws.state     = state;
  ws.params    = params;
  ws.gauge     = gauge->data;
  ws.tmp_e     = temps;
  ws.tmp2_e    = qx(step_even)(state, temps);
  ws.tmp_o     = qx(step_even)(state, temps);
  ws.flops     = &flops;
  ws.sent      = &sent;
  ws.received  = &received;

  lv = q(latvec_c_view)(deflator->dim, deflator->Ls, hfermion->even);
  q(df_inject)(deflator, &ws, lv);

  END_TIMING(state, flops, sent, received);
  if (ptr)
    q(free)(state, ptr, ptr_size);
  return 0;
}
