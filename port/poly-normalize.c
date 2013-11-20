#include <mdwf.h>
#include <math.h>

int 
Q(poly_normalize)(
            int poly_n, 
            double poly_a[],
            double poly_b[],
            double poly_c[],
            double x0,
            double tol)
{
    double p0, s0,
           p1, s1,
           p2, s2,
           aux;
    int i;
    if (poly_n <= 0 
            || 0 == poly_a
            || 0 == poly_b
            || 0 == poly_c)
        return 1;
    
    p0 = poly_c[0];
    s0 = (fabs(p0) < tol) ? 1. : 1. / p0;
    p1 = poly_a[0] + poly_b[0] * x0;
    s1 = (fabs(p1) < tol) ? 1. : 1. / p1;

    poly_c[0] *= s0;
    poly_a[0] *= s1;
    poly_b[0] *= s1;

    for (i = 1 ; i < poly_n ; i++) {
        p2 = (poly_a[i] + poly_b[i] * x0) * p1 + poly_c[i] * p0;
        s2 = (fabs(p2) < tol) ? 1. : 1. / p2;

        poly_a[i] *= (s2 / s1);
        poly_b[i] *= (s2 / s1);
        poly_c[i] *= (s2 / s0);

        p0 = p1;
        p1 = p2;

        s0 = s1;
        s1 = s2;
    }
    return 0;
}

