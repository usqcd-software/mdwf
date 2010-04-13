#include <stdio.h>
#include <stdlib.h>
#include <mdwf.h>

int
main(int argc, char *argv[])
{
    int size;
    int ls;
    int width;
    int f_size;
    int v_size;

    if (argc != 4) {
        fprintf(stderr, "Usage: check-vfsize size ls width\n");
        return 1;
    }
    size = atoi(argv[1]);
    ls = atoi(argv[2]);
    width = atoi(argv[3]);

    f_size = qx(sizeof_fermion)(size, ls);
    v_size = qx(sizeof_vfermion)(size, ls, width);

    printf("%d %d %d %d %d  %d %d\n",
           size, ls, width, v_size, f_size, v_size / f_size, v_size % f_size);



    return 0;
}
