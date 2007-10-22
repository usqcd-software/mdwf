#include <complex.h>
void
foo1(_Complex REAL *result,
     const _Complex REAL *a,
     const _Complex REAL *b)
{
    *result = *a;
}

void
foo2(_Complex REAL *result,
     const _Complex REAL *a,
     const _Complex REAL *b)
{
    *result = *a + *b;
}

void
foo3(_Complex REAL *result,
     const _Complex REAL *a,
     const _Complex REAL *b)
{
    *result = *a * *b;
}

void
foo4(_Complex REAL *result,
     const _Complex REAL *a,
     const _Complex REAL *b)
{
    *result = *a + I * *b;
}

void
foo5(_Complex REAL *result,
     const _Complex REAL *a,
     double s,
     const _Complex REAL *b)
{
    *result = *a + s * *b;
}
