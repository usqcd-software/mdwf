#include <qop-mdwf3.h>
#include "optest.h"

char *op_name = "Dirac Domain Wall operator";


void
setup_bc(void)
{
    int i;

    for (i = 0; i < lattice[4]; i++) {
	b5[i] = 1.0 / (i + 2.0);
	c5[i] = 1.0 / (lattice[4] + 1.0 - i * i);
    }
}

double
read_gauge(int dir,
	   const int pos[4],
	   int a, int b,
	   int re_im,
	   void *env)
{
    if (a == b && re_im == 0)
	return 1.0;
    else
	return 0.0;
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
	xprint("   fermion[%d,%d,%d,%d;%d][%d,%d].%d = %g",
	       pos[0], pos[1], pos[2], pos[3], pos[4],
	       c, d, re_im, value);
}

int operator(void)
{
    return QOP_MDWF_DDW_operator(result, params, gauge, fermion);
}