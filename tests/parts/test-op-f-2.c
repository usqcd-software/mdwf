#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "optest.h"
#include "op-routines.h"

char *op_name = "Operator lib:F";

double
read_gauge(int dir,
	   const int pos[4],
	   int a, int b,
	   int re_im,
	   void *env)
{
    int i;
    unsigned int s;
    double v;

    s = sum_init(seed_u);
    s = sum_add(s, dir);
    s = sum_add(s, a);
    for (i = 0; i < 4; i++)
	s = sum_add(s, pos[i]);
    s = sum_add(s, b);
    s = sum_add(s, re_im);
    s = sum_add(s, seed_u);
    v = sum_fini(s);

    xprint(" -u/2[%2d %2d %2d %2d].%1d(%1d,%1d)%d = %12.8f",
	   pos[0], pos[1], pos[2], pos[3], dir, a, b, re_im, -v/2);
    return v;
}

double
read_fermion(const int pos[5],
	     int c, int d,
	     int re_im,
	     void *env)
{
    int i;
    for (i = 0; i < 5; i++) {
	if (pos[i] != fermion_pos[i])
	    return 0.0;
    }
    if (c != fermion_color ||
	d != fermion_dirac ||
	re_im != fermion_reim)
	return 0.0;
    return 1.0;
}

void
write_fermion(const int pos[5],
	      int c, int d,
	      int re_im,
	      double value,
	      void *env)
{
    if (value != 0.0)
	xprint("   fermion[%d,%d,%d,%d;%d][%d,%d].%d = %12.8f",
	       pos[0], pos[1], pos[2], pos[3], pos[4],
	       c, d, re_im, value);
}

int operator(void)
{
    return op_F(result, gauge, fermion);
}
