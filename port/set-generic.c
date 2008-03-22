#include <mdwf.h>

int
Q(set_generic)(struct Q(Parameters) **params_ptr,
	       struct Q(State) *state,
	       const double b_5[],
	       const double c_5[],
	       double M_5,
	       double m)
{
  char *msg = NULL;
  int i;
  int Ls;
  struct Q(Parameters) *params;
  int abs;
  int iabs;
  double u, v, w;

#define CHECK(p, m) do { if ((p) == NULL) { msg = (m); goto error; }} while(0)

  if (state == NULL || state->error_latched)
    return 1;

  Ls = state->Ls;
  if (params_ptr == NULL)
    return q(set_error)(state, 0, "set_generic(): NULL pointer");

  *params_ptr = NULL;
  params = q(malloc)(state, sizeof (struct Q(Parameters)));
  CHECK(params, "set_generic(): Not enough memory for parameters");
  memset(params, 0, sizeof (struct Q(Parameters)));

  abs = q(sizeof_ABTable)(Ls);
  params->ATable = q(malloc)(state, abs);
  CHECK(params->ATable, "set_generic(): Not enough memory for A table");
#define u_a(j) (c_5[j] * (M_5 + 4) - 1)
#define v_a(j) (-m * u_a(j))
#define w_a(j) (b_5[j] * (M_5 + 4) + 1)

  u = u_a(0);
  v = v_a(0);
  w = w_a(0);
  q(put_ABTable)(params->ATable, 0, Ls-1, 1, w, v, u);
  for (i = 1; i < Ls - 1; i++) {
      u = u_a(i);
      w = w_a(i);
      q(put_ABTable)(params->ATable, i, i-1, i+1, w, u, u);
  }
  u = u_a(Ls-1);
  v = v_a(Ls-1);
  w = w_a(Ls-1);
  q(put_ABTable)(params->ATable, Ls-1, Ls-2, 0, w, u, v);

  params->AxTable = q(malloc)(state, abs);
  CHECK(params->AxTable, "set_generic(): Not enough memory for A* table");
  u = u_a(1);
  v = v_a(Ls-1);
  w = w_a(0);
  q(put_ABTable)(params->AxTable, 0, 1, Ls-1, w, u, v);
  for (i = 1; i < Ls - 1; i++) {
      double up = u_a(i+1);
      double um = u_a(i-1);
      w = w_a(i);
      q(put_ABTable)(params->AxTable, i, i+1, i-1, w, up, um);
  }
  u = u_a(Ls-2);
  v = v_a(0);
  w = w_a(Ls-1);
  q(put_ABTable)(params->AxTable, Ls-1, 0, Ls-2, w, v, u);
#undef u_a
#undef v_a
#undef w_a

  params->BTable = q(malloc)(state, abs);
  CHECK(params->BTable, "set_generic(): Not enough memory for B table");

#define u_b(j) (c_5[j])
#define v_b(j) (-m * u_b(j))
#define w_b(j) (b_5[j])

  u = u_b(0);
  v = v_b(0);
  w = w_b(0);
  q(put_ABTable)(params->BTable, 0, Ls-1, 1, w, v, u);
  for (i = 1; i < Ls - 1; i++) {
      u = u_b(i);
      w = w_b(i);
      q(put_ABTable)(params->BTable, i, i-1, i+1, w, u, u);
  }
  u = u_b(Ls-1);
  v = v_b(Ls-1);
  w = w_b(Ls-1);
  q(put_ABTable)(params->BTable, Ls-1, Ls-2, 0, w, u, v);

  params->BxTable = q(malloc)(state, abs);
  CHECK(params->BxTable, "set_generic(): Not enough memory for B* table");
  u = u_b(1);
  v = v_b(Ls-1);
  w = w_b(0);
  q(put_ABTable)(params->BxTable, 0, 1, Ls-1, w, u, v);
  for (i = 1; i < Ls - 1; i++) {
      double up = u_b(i+1);
      double um = u_b(i-1);
      w = w_b(i);
      q(put_ABTable)(params->BxTable, i, i+1, i-1, w, up, um);
  }
  u = u_b(Ls-2);
  v = v_b(0);
  w = w_b(Ls-1);
  q(put_ABTable)(params->BxTable, Ls-1, 0, Ls-2, w, v, u);
#undef u_b
#undef v_b
#undef w_b

  /* XXX fill iATable and iBTable */
  iabs = q(sizeof_ABiTable)(Ls);
  params->AipTable = q(malloc)(state, iabs);
  CHECK(params->AipTable, "set_generic(): Not enough memory for 1/A + table");
  params->AimTable = q(malloc)(state, iabs);
  CHECK(params->AimTable, "set_generic(): Not enough memory for 1/A - table");
  params->BipTable = q(malloc)(state, iabs);
  CHECK(params->BipTable, "set_generic(): Not enough memory for 1/B + table");
  params->BimTable = q(malloc)(state, iabs);
  CHECK(params->BimTable, "set_generic(): Not enough memory for 1/B - table");
  q(put_ABiTableZ)(params->AipTable, 0.0);
  q(put_ABiTableZ)(params->AimTable, 0.0);
  q(put_ABiTableZ)(params->BipTable, 0.0);
  q(put_ABiTableZ)(params->BimTable, 0.0);
  for (i = 1; i < Ls; i++) {
      q(put_ABiTable)(params->AipTable, i,  0.0,  0.0, 0.0);
      q(put_ABiTable)(params->AimTable, i,  0.0,  0.0, 0.0);
      q(put_ABiTable)(params->BipTable, i,  0.0,  0.0, 0.0);
      q(put_ABiTable)(params->BimTable, i,  0.0,  0.0, 0.0);
  }  
  /* XXX */

  BEGIN_TIMING(state);
  params->state = state;
  *params_ptr = params;
  END_TIMING(state, 0, 0, 0);

  return 0;
#undef CHECK
 error:
  Q(free_parameters)(&params);
  return q(set_error)(state, 0, msg);
}
