#include <mdwf.h>

void *
q(malloc)(struct Q(State) *state, size_t bytes)
{
  void *ptr = malloc(bytes);

  if (ptr == NULL)
    return NULL;

  if (state) {
    state->allocated += bytes;
    if (state->allocated > state->max_allocated) {
      state->max_allocated = state->allocated;
    }
  }
  return ptr;
}
