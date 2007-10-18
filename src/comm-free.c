#include <mdwf.h>

static void
eo_free(struct eo_lattice *eo, struct Q(State) *state)
{
  int k;

  if (eo->real_size == 0)
    return;

  for (k = eo->hrt_count; k--;)
    QMP_free_msghandle(eo->hrt[k]);
  eo->hrt_count = 0;

  for (k = eo->hsv_count; k--;)
    QMP_free_msghandle(eo->hsv[k]);
  eo->hsv_count = 0;

  for (k = eo->mh_count; k--;)
    QMP_free_msgmem(eo->mh[k]);
  eo->mh_count = 0;

  for (k = eo->mem_count; k--;)
    QMP_free_memory(eo->mem[k]);
  eo->mem_count = 0;

  for (k = 0; k < Q(DIM); k++) {
    eo->send_up_buf[k] = 0;
    eo->send_down_buf[k] = 0;
    eo->receive_buf[k] = 0;
    eo->receive_buf[k + Q(DIM)] = 0;
  }
  eo->real_size = 0;
}

int
q(free_comm)(struct Q(State) *state)
{
  eo_free(&state->even, state);
  eo_free(&state->odd, state);
  return 1;
}
