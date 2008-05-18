#include <mdwf.h>

const char *
Q(error)(struct Q(State) *state)
{
  const char *msg;

  if (state == 0)
    return "NULL MDWF state";

  msg = state->error;
  if (state->fatal_error == 0) {
    state->error_latched = 0;
    state->error = NULL;
  }
  return msg;
}

int
q(set_error)(struct Q(State) *state, int fatal, const char *error)
{
    if (state->error_latched == 0) {
	state->error_latched = 1;
	state->error = error;
	state->fatal_error = fatal;
    }
    return 1;
}
