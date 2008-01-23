#include <mdwf.h>

static void
free_atable(struct Q(State) *state, struct ABTable **table)
{
    int Ls = state->Ls;
    int abs = q(sizeof_ABTable(Ls));

    if (*table == 0)
	return;
    q(free)(state, *table, abs);
    *table = 0;
}

static void
free_aitable(struct Q(State) *state, struct ABiTable **table)
{
    int Ls = state->Ls;
    int abis = q(sizeof_ABiTable(Ls));

    if (*table == 0)
	return;
    q(free)(state, *table, abis);
    *table = 0;
}

void
Q(free_parameters)(struct Q(Parameters) **parameters_ptr)
{
  struct Q(State) *state;

  if (parameters_ptr == 0 || *parameters_ptr == 0)
    return;

  state = (*parameters_ptr)->state;
  BEGIN_TIMING(state);
  free_atable(state, &(*parameters_ptr)->ATable);
  free_atable(state, &(*parameters_ptr)->BTable);
  free_aitable(state, &(*parameters_ptr)->AipTable);
  free_aitable(state, &(*parameters_ptr)->AimTable);
  free_aitable(state, &(*parameters_ptr)->BipTable);
  free_aitable(state, &(*parameters_ptr)->BimTable);
  END_TIMING(state, 0, 0, 0);
  q(free)(state, *parameters_ptr, sizeof (struct Q(Parameters)));
  *parameters_ptr = 0;
}
