#include <mdwf.h>

void
q(free)(struct Q(State) *state, void *ptr, size_t bytes)
{
    if (ptr == NULL)
	return;
    
    free(ptr);
    
    if (state) {
	state->allocated -= bytes;
	Q(fini)(&state);
    }
}
