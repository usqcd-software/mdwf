#include <mdwf.h>

int
Q(performance)(double *time_sec,
	       long long *flops,
	       long long *sent,
	       long long *received,
	       struct Q(State) *state)
{
  if (state == NULL || state->error_latched)
    return 1;

  if (time_sec)
    *time_sec = state->time_sec;
  if (flops)
    *flops = state->flops;
  if (sent)
    *sent = state->sent;
  if (received)
    *received = state->received;
  return 0;
}
