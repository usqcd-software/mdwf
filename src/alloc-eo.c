#include <mdwf.h>

void *
q(allocate_eo)(struct Q(State) *state,
               size_t *size, void **aligned_ptr,
               size_t hdr_size, int even_count, int odd_count, size_t fsize)
{
#define HFSIZE(n,s4) ((n) * (state->Ls * fsize * (s4) + CACHE_LINE_SIZE - 1))
  size_t total_size = hdr_size
                    + HFSIZE(even_count, state->even.full_size)
                    + HFSIZE(odd_count, state->odd.full_size);

  return q(allocate_aligned)(state, size, aligned_ptr, hdr_size, total_size);
}
