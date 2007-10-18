#include <mdwf.h>

void *
q(allocate_eo)(struct Q(State) *state,
               size_t *size, void **aligned_ptr,
               size_t hdr_size, int even_count, int odd_count, size_t fsize)
{
  int es;
  int os;
  size_t total_size;

  qx(sizeof_fermion)(&es, state->even.full_size, state->Ls);
  qx(sizeof_fermion)(&os, state->odd.full_size, state->Ls);

#define HFSIZE(n,xx) ((n) * ((xx) + CACHE_LINE_SIZE - 1))
  total_size = hdr_size + HFSIZE(even_count, es) + HFSIZE(odd_count, os);

  return q(allocate_aligned)(state, size, aligned_ptr, hdr_size, total_size);
}
