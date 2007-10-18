#include <mdwf.h>

void *
q(allocate_aligned)(struct Q(State) *state,
                    size_t *size, void **aligned_ptr,
                    size_t hdr_size, size_t bulk_size)
{
  size_t total_size = hdr_size + bulk_size + CACHE_LINE_SIZE - 1;
  void *ptr;

  if (state == NULL || size == NULL || aligned_ptr == NULL)
    return NULL;
  ptr = q(malloc)(state, total_size);
  if (ptr == NULL) {
    *size = 0;
    *aligned_ptr = NULL;
  } else {
    *size = total_size;
    *aligned_ptr = ALIGN(ptr + hdr_size);
  }
  return ptr;
}
