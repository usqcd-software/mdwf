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
  free_atable(state, &(*parameters_ptr)->ApTable);
  free_atable(state, &(*parameters_ptr)->AmTable);
  free_atable(state, &(*parameters_ptr)->AxpTable);
  free_atable(state, &(*parameters_ptr)->AxmTable);
  free_atable(state, &(*parameters_ptr)->BpTable);
  free_atable(state, &(*parameters_ptr)->BmTable);
  free_atable(state, &(*parameters_ptr)->BxpTable);
  free_atable(state, &(*parameters_ptr)->BxmTable);
  free_aitable(state, &(*parameters_ptr)->AipTable);
  free_aitable(state, &(*parameters_ptr)->AimTable);
  free_aitable(state, &(*parameters_ptr)->BipTable);
  free_aitable(state, &(*parameters_ptr)->BimTable);
  free_aitable(state, &(*parameters_ptr)->AxipTable);
  free_aitable(state, &(*parameters_ptr)->AximTable);
  free_aitable(state, &(*parameters_ptr)->BxipTable);
  free_aitable(state, &(*parameters_ptr)->BximTable);
  END_TIMING(state, 0, 0, 0);
  q(free)(state, *parameters_ptr, sizeof (struct Q(Parameters)));
  *parameters_ptr = 0;
}
