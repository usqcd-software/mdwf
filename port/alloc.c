#include <mdwf.h>

void *
q(malloc)(struct Q(State) *state, size_t bytes)
{
    void *ptr = malloc(bytes);
    
    if (ptr == NULL)
	return NULL;
    
    if (state) {
	state->used++;
	state->allocated += bytes;
	if (state->allocated > state->max_allocated) {
	    state->max_allocated = state->allocated;
	}
	printf("%s(%d): %p, state=%p, used=%d, saved=%d\n",
	       __FUNCTION__, (int)bytes,
	       ptr,
	       state,
	       state->used, state->saved);
    } else {
	printf("%s(%d): %p, NULL\n",
	       __FUNCTION__, (int)bytes,
	       ptr);
    }
    return ptr;
}
