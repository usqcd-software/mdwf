#include <mdwf.h>

int
Q(performance)(double *time_sec,
	       long long *flops,
	       long long *sent,
	       long long *received,
	       struct Q(State) *state)
{
  if (state == 0 || state->error_latched)
    return 1;

  if (time_sec == 0 || flops == 0 || sent == 0 || received == 0) {
    q(set_error)(state, 0, "performance(): NULL pointer(s)");
    return 1;
  }
  
  *time_sec = state->time_sec;
  *flops = state->flops;
  *sent = state->sent;
  *received = state->received;
  return 0;
}
