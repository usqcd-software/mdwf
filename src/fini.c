#include <mdwf.h>

static void
free_eo(struct Q(State) *state, struct eo_lattice *eo)
{
  int i;
  int ns;
  int us;
  int ds;

  for (i = 0; i < Q(DIM); i++) {
    q(sizeof_up_pack)(&us, eo->send_up_size[i]);
    if (eo->up_pack[i])
      q(free)(state, eo->up_pack[i], us);
    q(sizeof_down_pack)(&ds, eo->send_down_size[i]);
    if (eo->down_pack[i])
      q(free)(state, eo->down_pack[i], ds);
  }

  q(sizeof_neighbor)(&ns, eo->full_size);
  if (eo->body_neighbor)
    q(free)(state, eo->body_neighbor, ns);

  if (eo->lx2v)
    q(free)(state, eo->lx2v, eo->full_size * sizeof (int));

  if (eo->v2lx)
    q(free)(state, eo->v2lx, state->volume * sizeof (int));

  memset(eo, 0, sizeof (struct eo_lattice));
}

void
q(cleanup_state)(struct Q(State) *state)
{
  q(free_comm)(state);
  free_eo(state, &state->even);
  free_eo(state, &state->odd);
  if (state->lx2v)
    q(free)(state, state->lx2v, state->volume * sizeof (int));
  if (state->v2lx)
    q(free)(state, state->v2lx, state->volume * sizeof (int));
  memset(state, 0, sizeof (struct Q(State)));
}

void
Q(fini)(struct Q(State) **state_ptr)
{
  struct Q(State) *state;

  if (state_ptr == 0 || *state_ptr == 0)
    return;

  state = *state_ptr;
  *state_ptr = NULL;
  state->used--;
  if (state->used != 0)
    return;

  q(cleanup_state)(state);
  free(state);
}
