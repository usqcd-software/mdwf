/* Browsing the environment */

#include <stdlib.h>
#include <stdio.h>

extern char **environ;
int
main(int argc, char *argv[])
{
    char **e;
    for (e = environ; *e; e++)
	printf("[%p] %s\n", e, *e);
    return 0;
}
