#include <mdwf.h>

static void
free_eo(struct Q(State) *state, struct eo_lattice *eo)
{
  int i;

  if (eo->layout2vector)
    free(eo->layout2vector);
  if (eo->body_neighbor)
    free(eo->body_neighbor);
  /* .face_neighbor was allocaed with .body_neighbor */

  for (i = 0; i < Q(DIM); i++) {
    if (eo->up_pack[i])
      free(eo->up_pack[i]);
    if (eo->down_pack[i])
      free(eo->down_pack[i]);
  }
  memset(eo, 0, sizeof (struct eo_lattice));
}

void
Q(fini)(struct Q(State) **state_ptr)
{
  if (state_ptr == 0 || *state_ptr == 0)
    return;

  if (--((*state_ptr)->used))
    return;

  q(free_comm)(*state_ptr);
  free_eo(*state_ptr, &(*state_ptr)->even);
  free_eo(*state_ptr, &(*state_ptr)->odd);
  q(free)(*state_ptr, (*state_ptr)->layout2linear,
	  (*state_ptr)->volume * sizeof (int));
  free(*state_ptr);
  *state_ptr = 0;
}
