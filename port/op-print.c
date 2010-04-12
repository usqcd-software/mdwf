#include <mdwf.h>
#include <stdarg.h>

void
qx(zprint)(struct Q(State) *state,
           const char *source,
           const char *fmt,
           ...)
{
    if (state->master_p) {
        char buffer[4096];
        va_list va;
        int len;

        snprintf(buffer, sizeof (buffer) - 1, "MDWF: %s(%c): ",
                 source, Q(DEFAULT_PRECISION));
        len = strlen(buffer);
        va_start(va, fmt);
        vsnprintf(buffer + len, sizeof (buffer) - len - 1, fmt, va);
        va_end(va);
        printf("%s\n", buffer);
    }
}
