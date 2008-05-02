#include <mdwf.h>
#include <stdarg.h>

void
qx(zprint)(struct Q(State) *state,
	   const char *source,
	   const char *fmt,
	   ...)
{
    char buffer[4096];
    va_list va;

    if (state->master_p) {
	va_start(va, fmt);
	vsnprintf(buffer, sizeof (buffer) - 1, fmt, va);
	va_end(va);
	QMP_printf("MDWF: %s(%c): %s\n", source, Q(DEFAULT_PRECISION), buffer);
    }
}
