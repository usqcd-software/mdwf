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
  u = c_5[0]*(M_5 + 4) - 1;
  v = -m * u;
  w = b_5[0]*(M_5 + 4) - 1;
/* XXX q(put_ABTable)(params->ATable, 0, Ls-1, 1, w, v, u); */
  q(put_ABTable)(params->ATable, 0, Ls-1, 1, 0.0, 0.0, 0.0);
  for (i = 1; i < Ls - 1; i++) {
    u = c_5[i]*(M_5 + 4) - 1;
    w = b_5[i]*(M_5 + 4) - 1;
/* XXX    q(put_ABTable)(params->ATable, i, i-1, i+1, w, u, u); */
    q(put_ABTable)(params->ATable, i, i-1, i+1, 0.0, 0.0, 0.0);
  }
  u = c_5[Ls-1]*(M_5 + 4) - 1;
  v = -m * u;
  w = b_5[Ls-1]*(M_5 + 4) - 1;
/* XXX  q(put_ABTable)(params->ATable, Ls-1, Ls-2, 0, w, u, v); */
  q(put_ABTable)(params->ATable, Ls-1, Ls-2, 0, 0.0, 0.0, 0.0);

  params->BTable = q(malloc)(state, abs);
  CHECK(params->BTable, "set_generic(): Not enough memory for B table");
  u = c_5[0];
  v = -m * u;
  w = b_5[0];
/* XXX  q(put_ABTable)(params->BTable, 0, Ls-1, 1, w, v, u); */
  q(put_ABTable)(params->BTable, 0, Ls-1, 1, 0.0, 0.0, 0.0);
  for (i = 1; i < Ls - 1; i++) {
    u = c_5[i];
    w = b_5[i];
/* XXX   q(put_ABTable)(params->BTable, i, i-1, i+1, w, u, u); */
    q(put_ABTable)(params->BTable, i, i-1, i+1, 0.0, 0.0, 0.0);
  }
  u = c_5[Ls-1];
  v = -m * u;
  w = b_5[Ls-1];
/* XXX  q(put_ABTable)(params->BTable, Ls-1, Ls-2, 0, w, u, v); */
  q(put_ABTable)(params->BTable, Ls-1, Ls-2, 0, 0.0, 0.0, 0.0);

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
