#include <complex.h>

void
foo(int n, void *va, void *vb)
{
   int i;
   double _Complex *a = (double _Complex *)va;
   double _Complex *b = (double _Complex *)vb;
   double _Complex c;

#pragma disjoint(*a,*b)
   __alignx(16, a);
   __alignx(16, b);
   c = *a;
   for (i = 0; i < n; i++)
     c = *a++ + c * *b++;
   *a = c;
}

