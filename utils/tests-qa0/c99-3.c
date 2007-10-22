#include <complex.h>

/* a = b */
void
foo1(double _Complex *a,
     const double _Complex *b)
{
    *a = *b;
}

/* a = -b */
void
foo2(double _Complex *a,
     const double _Complex *b)
{
    *a = -*b;
}

/* a = I * b */
void
foo3(double _Complex *a,
     const double _Complex *b)
{
    *a = I * *b;
}

/* a = -I * b */
void
foo4(double _Complex *a,
     const double _Complex *b)
{
    *a = - I * *b;
}

/* a = b + c */
void
foo5(double _Complex *a,
     const double _Complex *b,
     const double _Complex *c)
{
    *a = *b + *c;
}

/* a = b - c */
void
foo6(double _Complex *a,
     const double _Complex *b,
     const double _Complex *c)
{
    *a = *b - *c;
}

/* a = b * c */
void
foo7(double _Complex *a,
     const double _Complex *b,
     const double _Complex *c)
{
    *a = *b * *c;
}

/* a = conj(b) * c */
void
foo8(double _Complex *a,
     const double _Complex *b,
     const double _Complex *c)
{
    *a = conj(*b) * *c;
}

/* a = d + b * c */
void
foo9(double _Complex *a,
     const double _Complex *d,
     const double _Complex *b,
     const double _Complex *c)
{
    *a = *d + *b * *c;
}

/* a = d + conj(b) * c */
void
foo10(double _Complex *a,
      const double _Complex *d,
      const double _Complex *b,
      const double _Complex *c)
{
    *a = *d + conj(*b) * *c;
}

/* a = d + I * c */
void
foo11(double _Complex *a,
     const double _Complex *d,
     const double _Complex *c)
{
    *a = *d - I * *c;
}

/* a = d - I * c */
void
foo12(double _Complex *a,
      const double _Complex *d,
      const double _Complex *c)
{
    *a = *d - I * *c;
}

/* a = d + s * c */
void
foo13(double _Complex *a,
      const double _Complex *d,
      double s,
      const double _Complex *c)
{
    *a = *d + s * *c;
}
