#include <mdwf.h>

static void
free_eo(struct Q(State) *state, struct eo_lattice *eo)
{
  int i;
  int up_pack_size;
  int neighbor_size;

  q(sizeof_up_pack)(&up_pack_size);
  q(sizeof_neighbor)(&neighbor_size);
  for (i = 0; i < Q(DIM); i++) {
    if (eo->up_pack[i])
      q(free)(state, eo->up_pack[i], eo->send_up_size[i] * up_pack_size);
    if (eo->down_pack[i])
      q(free)(state, eo->down_pack[i], eo->send_down_size[i] * sizeof (int));
  }
  if (eo->body_neighbor)
    q(free)(state, eo->body_neighbor, eo->body_size * neighbor_size);

  if (eo->layout2vector)
    q(free)(state, eo->layout2vector, eo->full_size * sizeof (int));

  memset(eo, 0, sizeof (struct eo_lattice));
}

void
q(cleanup_state)(struct Q(State) *state)
{
  q(free_comm)(state);
  free_eo(state, &state->even);
  free_eo(state, &state->odd);
  q(free)(state, state->layout2vector,
	  state->volume * sizeof (int));
  memset(state, 0, sizeof (struct Q(State)));
}

void
Q(fini)(struct Q(State) **state_ptr)
{
  if (state_ptr == 0 || *state_ptr == 0)
    return;

  if (--((*state_ptr)->used))
    return;

  q(cleanup_state)(*state_ptr);
  free(*state_ptr);
  *state_ptr = 0;
}
