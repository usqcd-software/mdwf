#include <mdwf.h>

void
q(free)(struct Q(State) *state, void *ptr, size_t bytes)
{
    if (ptr == NULL)
	return;
    
    free(ptr);
    
    if (state) {
	printf("%s(%p,%d): state=%p, used=%d, saved=%d\n",
	       __FUNCTION__, ptr, (int)bytes,
	       state,
	       state->used, state->saved);
	state->allocated -= bytes;
	Q(fini)(&state);
    }
}
