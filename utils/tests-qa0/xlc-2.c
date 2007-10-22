#include <complex.h>

/* a = b */
void
foo1(double _Complex *a,
     const double _Complex *b)
{
#pragma disjoint(*a,*b)
    __alignx(16, a);
    __alignx(16, b);
    *a = *b;
}

/* a = -b */
void
foo2(double _Complex *a,
     const double _Complex *b)
{
#pragma disjoint(*a,*b)
    __alignx(16, a);
    __alignx(16, b);
    *a = __fpneg(*b);
}

/* a = I * b */
void
foo3(double _Complex *a,
     const double _Complex *b)
{
#pragma disjoint(*a,*b)
    __alignx(16, a);
    __alignx(16, b);
    double _Complex j;

    j = __cmplx(-1.0, 1.0);
    *a = __fxmul(*b, j);
}

/* a = -I * b */
void
foo4(double _Complex *a,
     const double _Complex *b)
{
#pragma disjoint(*a,*b)
    __alignx(16, a);
    __alignx(16, b);
    double _Complex j;

    j = __cmplx(1.0, -1.0);
    *a = __fxmul(*b,j);
}

/* a = b + c */
void
foo5(double _Complex *a,
     const double _Complex *b,
     const double _Complex *c)
{
#pragma disjoint(*a,*b)
#pragma disjoint(*a,*c)
    __alignx(16, a);
    __alignx(16, b);
    __alignx(16, c);

    *a = __fpadd(*b,*c);
}

/* a = b - c */
void
foo6(double _Complex *a,
     const double _Complex *b,
     const double _Complex *c)
{
#pragma disjoint(*a,*b)
#pragma disjoint(*a,*c)
    __alignx(16, a);
    __alignx(16, b);
    __alignx(16, c);

    *a = __fpsub(*b,*c);
}

/* a = b * c */
void
foo7(double _Complex *a,
     const double _Complex *b,
     const double _Complex *c)
{
    double _Complex bc;
#pragma disjoint(*a,*b)
#pragma disjoint(*a,*c)
    __alignx(16, a);
    __alignx(16, b);
    __alignx(16, c);

    bc = __fxpmul(*c, __creal(*b));
    *a = __fxcpnpma(bc, *c, __cimag(*b));
}

/* a = conj(b) * c */
void
foo8(double _Complex *a,
     const double _Complex *b,
     const double _Complex *c)
{
    double _Complex bc;
#pragma disjoint(*a,*b)
#pragma disjoint(*a,*c)
    __alignx(16, a);
    __alignx(16, b);
    __alignx(16, c);

    bc = __fxpmul(*c, __creal(*b));
    *a = __fxcpnsma(bc, *c, __cimag(*b));
}

/* a = d + b * c */
void
foo9(double _Complex *a,
     const double _Complex *d,
     const double _Complex *b,
     const double _Complex *c)
{
    double _Complex bc;
#pragma disjoint(*a,*b)
#pragma disjoint(*a,*c)
#pragma disjoint(*a,*d)
    __alignx(16, a);
    __alignx(16, b);
    __alignx(16, c);
    __alignx(16, d);

    bc = __fxcpmadd(*d,*c, __creal(*b));
    *a = __fxcpnpma(bc, *c, __cimag(*b));
}

/* a = d + conj(b) * c */
void
foo10(double _Complex *a,
      const double _Complex *d,
      const double _Complex *b,
      const double _Complex *c)
{
    double _Complex bc;
#pragma disjoint(*a,*b)
#pragma disjoint(*a,*c)
#pragma disjoint(*a,*d)
    __alignx(16, a);
    __alignx(16, b);
    __alignx(16, c);
    __alignx(16, d);

    bc = __fxcpmadd(*d, *c, __creal(*b));
    *a = __fxcpnsma(bc, *c, __cimag(*b));
}

/* a = d + I * c */
void
foo11(double _Complex *a,
     const double _Complex *d,
     const double _Complex *c)
{
#pragma disjoint(*a,*c)
#pragma disjoint(*a,*d)
    __alignx(16, a);
    __alignx(16, c);
    __alignx(16, d);

    *a = __fxcpnpma(*d, *c, 1.0);
}

/* a = d - I * c */
void
foo12(double _Complex *a,
      const double _Complex *d,
      const double _Complex *c)
{
#pragma disjoint(*a,*c)
#pragma disjoint(*a,*d)
    __alignx(16, a);
    __alignx(16, c);
    __alignx(16, d);

    *a = __fxcpnsma(*d, *c, 1.0);
}

/* a = d + s * c */
void
foo13(double _Complex *a,
      const double _Complex *d,
      double s,
      const double _Complex *c)
{
#pragma disjoint(*a,*c)
#pragma disjoint(*a,*d)
    __alignx(16, a);
    __alignx(16, c);
    __alignx(16, d);

    *a = __fxcpnpma(*d, *c, s);
}
