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

  abs = q(sizeof_ABTable)(state->Ls);
  params->ATable = q(malloc)(state, abs);
  CHECK(params->ATable, "set_generic(): Not enough memory for A table");
  u = c_5[0]*(M_5 + 4) - 1;
  v = -m * u;
  w = b_5[0]*(M_5 + 4) - 1;
  q(put_ABTable)(params->ATable, 0, Ls-1, 1, w, v, u);
  for (i = 1; i < Ls - 1; i++) {
    u = c_5[i]*(M_5 + 4) - 1;
    w = b_5[i]*(M_5 + 4) - 1;
    q(put_ABTable)(params->ATable, i, i-1, i+1, w, u, u);
  }
  u = c_5[Ls-1]*(M_5 + 4) - 1;
  v = -m * u;
  w = b_5[Ls-1]*(M_5 + 4) - 1;
  q(put_ABTable)(params->ATable, Ls-1, Ls-2, 0, w, u, v);

  params->BTable = q(malloc)(state, abs);
  CHECK(params->BTable, "set_generic(): Not enough memory for B table");
  u = c_5[0];
  v = -m * u;
  w = b_5[0];
  q(put_ABTable)(params->BTable, 0, Ls-1, 1, w, v, u);
  for (i = 1; i < Ls - 1; i++) {
    u = c_5[i];
    w = b_5[i];
    q(put_ABTable)(params->BTable, i, i-1, i+1, w, u, u);
  }
  u = c_5[Ls-1];
  v = -m * u;
  w = b_5[Ls-1];
  q(put_ABTable)(params->BTable, Ls-1, Ls-2, 0, w, u, v);

  /* XXX */

  BEGIN_TIMING(state);
  state->used++;
  params->state = state;
  *params_ptr = params;
  END_TIMING(state, 0, 0, 0);

  return 0;
#undef CHECK
 error:
  if (params == NULL)
    return q(set_error)(state, 0, msg);

  if (params->ATable != NULL)
    q(free)(state, params->ATable, abs);
  if (params->BTable != NULL)
    q(free)(state, params->BTable, abs);
  q(free)(state, params, sizeof (struct Q(Parameters)));
  state->used--;
  return q(set_error)(state, 0, msg);
}
