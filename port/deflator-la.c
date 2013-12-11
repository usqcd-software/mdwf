#include <assert.h>
#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>
#include <qmp.h>

/* replace all references to `dim' and `Ls' with these macros */
#define DEFLATOR_DIM(pstate) ((pstate)->even.full_size)
#define DEFLATOR_LS(pstate) ((pstate)->Ls)

/* TODO replace float & double versions with dry code
    q(latvec_c_XXX), q(latvec_z_XXX) -> qx(latmat_XXX)
    q(latmat_c_XXX), q(latmat_z_XXX) -> qx(latmat_XXX)
 */

/* allocate & free */
latvec_c
q(latvec_c_alloc)(struct Q(State) *state)
{
    latvec_c res;
    res.dim     = DEFLATOR_DIM(state);
    res.Ls      = DEFLATOR_LS(state);
    res.mem_ptr = qx(allocate_eo)(state, &res.mem_size,
                                  (void *)&res.f, 0, 1, 0);
    if (res.mem_ptr == NULL)
        res.f = NULL;
    return res;
}
latvec_c 
q(latvec_c_view)(struct Q(State) *state, struct FermionF *f)
{
    latvec_c res;
    res.dim     = DEFLATOR_DIM(state);
    res.Ls      = DEFLATOR_LS(state);
    res.f       = f;
    res.mem_ptr = NULL;
    res.mem_size = 0;
    return res;
}
void 
q(latvec_c_copy)(latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_c_is_null(&x));
    assert(!latvec_c_is_null(&y));
    qx(f_copy)(y.f, y.dim, y.Ls, x.f);
}
void
q(latvec_c_zero)(latvec_c x)
{
    assert(!latvec_c_is_null(&x));
    qx(f_zero)(x.f, x.dim, x.Ls);
}
void
q(latvec_c_free)(struct Q(State) *state, latvec_c *v)
{
    if (NULL != v && NULL != v->mem_ptr && NULL != v->f) {
        q(free)(state, v->mem_ptr, v->mem_size);
        v->mem_ptr = NULL;
        v->f = NULL;
    }
}


doublecomplex 
q(lat_c_dotu)(latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_c_is_null(&x));
    assert(!latvec_c_is_null(&y));

    doublecomplex res = { 0., 0. };
    double s[2];

    qx(f_dot)(&s[0], &s[1], x.dim, x.Ls, x.f, y.f);
    QMP_sum_double_array(s, 2);
    res.r = s[0];
    res.i = s[1];
    
    return res;
}
void 
q(lat_c_scal_d)(double alpha, latvec_c x)
{
    assert(!latvec_c_is_null(&x));
    qx(f_rmul1)(x.f, x.dim, x.Ls, alpha);
}

void 
q(lat_c_axpy_d)(double alpha, latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_c_is_null(&x));
    assert(!latvec_c_is_null(&y));

    qx(f_add2)(y.f, y.dim, y.Ls, alpha, x.f);
}
double 
q(lat_c_nrm2)(latvec_c x) 
{
    assert(!latvec_c_is_null(&x));
    
    double res = 0.;
    qx(f_norm)(&res, x.dim, x.Ls, x.f);
    QMP_sum_double(&res);
    
    return res;
}

#if 0
latvec_z
q(latvec_z_alloc)(struct Q(State) *state, int dim)
{
    latvec_z res;
    res.dim     = dim;
    res.f       = NULL;
    /* TODO LATER allocate f */NOT_IMPLEMENTED;
    return res;
}
latvec_z
q(latvec_z_view)(int dim, struct FermionD *f)
{
    latvec_z res;
    res.dim     = dim;
    res.f       = f;
    return res;
}
void 
q(latvec_z_copy)(latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_z_is_null(&x));
    assert(!latvec_z_is_null(&y));
    /* TODO LATER copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}
void
q(latvec_z_free)(struct Q(State) *state, latvec_z *v)
{
    if (NULL != v && NULL != v->f) {
        /* TODO LATER free fermion */NOT_IMPLEMENTED;
        v->f = NULL;
    }
}

void 
q(latvec_zc_copy)(latvec_z x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_z_is_null(&x));
    assert(!latvec_c_is_null(&y));
    /* TODO LATER copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}
void 
q(latvec_cz_copy)(latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_c_is_null(&x));
    assert(!latvec_z_is_null(&y));
    /* TODO LATER copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}


doublecomplex 
q(lat_cz_dotu)(latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_c_is_null(&x));
    assert(!latvec_z_is_null(&y));

    doublecomplex res = { 0., 0. };
    /* TODO LATER return x^H . y */NOT_IMPLEMENTED;
    
    return res;
}
doublecomplex 
q(lat_z_dotu)(latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_z_is_null(&x));
    assert(!latvec_z_is_null(&y));

    doublecomplex res = { 0., 0. };
    /* TODO LATER return x^H . y */NOT_IMPLEMENTED;
    
    return res;
}
/* norm */
double 
q(lat_z_nrm2)(latvec_z x) 
{
    assert(!latvec_z_is_null(&x));
    double res = 0.;
    
    /* TODO LATER check that all is correct */NOT_IMPLEMENTED;
    qop_d3_clover_f_norm(&res, x.dim, x.f);
    res = res * res;
    
    return res;
}


void 
q(lat_c_scal)(doublecomplex alpha, latvec_c x)
{
    assert(!latvec_c_is_null(&x));
    /* TODO LATER implement x <- alpha * x */NOT_IMPLEMENTED;
}
void 
q(lat_c_axpy)(doublecomplex alpha, latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_c_is_null(&x));
    asseet(!latvec_c_is_null(&y));

    qx(f_add2)(y.f, y.dim, alpha, x.f);
}
void 
q(lat_cz_axpy)(doublecomplex alpha, latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_c_is_null(&x));
    asseert(!latvec_z_is_null(&y));
    /* TODO LATER implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
void 
q(lat_z_axpy)(doublecomplex alpha, latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_z_is_null(&x));
    assert(!latvec_z_is_null(&y));
    /* TODO LATER implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
void 
q(lat_cz_axpy_d)(double alpha, latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_c_is_null(&x));
    assert(!latvec_z_is_null(&y));
    /* TODO LATER implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
void 
q(lat_z_axpy_d)(double alpha, latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_z_is_null(&x));
    assert(!latvec_z_is_null(&y));
    /* TODO LATER implement y <- y + alpha * x */NOT_IMPLEMENTED;
}

#endif
    
    
latmat_c 
q(latmat_c_alloc)(struct Q(State) *state, int ncol)
{
    latmat_c res;
    res.dim     = DEFLATOR_DIM(state);
    res.Ls      = DEFLATOR_LS(state);
    res.size    = ncol;
    res.begin   = 0;
    res.len     = ncol;
    res.fv      = NULL;
    res.stride  = qx(strideof_vfermion)(res.dim, res.Ls);
    res.mem_ptr = q(allocate_aligned)(state, &res.mem_size, (void *)&res.fv, 0,
                                      qx(sizeof_vfermion)(res.dim, res.Ls, ncol));
    if (res.mem_ptr == NULL)
        res.fv = NULL;

    return res;
}

int
q(latmat_convert_to_blas)(latmat_c *m) 
{ 
    return 0; /* stub; modify if perform any transformations */
}
int
q(latmat_convert_from_blas)(latmat_c *m) 
{ 
    return 0; /* stub; modify if perform any transformations */
}

latmat_c 
q(latmat_c_view)(struct Q(State) *state, int size, struct vFermion *fv)
{
    latmat_c res;
    res.dim     = DEFLATOR_DIM(state);
    res.Ls      = DEFLATOR_LS(state);
    res.size    = size;
    res.begin   = 0;
    res.len     = size;
    res.fv       = fv;
    res.stride  = qx(strideof_vfermion)(res.dim, res.Ls);
    res.mem_ptr = NULL;
    res.mem_size = 0;
    return res;
}
void
q(latmat_c_free)(struct Q(State) *state, latmat_c *m) 
{
    if (NULL != m && NULL != m->mem_ptr && NULL != m->fv) {
        assert(0 == m->begin);
        assert(m->size == m->len);
        /* at least some chance that this is not a sub-matrix */
        void *mem_ptr = m->mem_ptr;
        size_t mem_size = m->mem_size;
        q(free)(state, mem_ptr, mem_size);
        if (m != mem_ptr) { 
            /* do not set default values if `m' is at mem_ptr */
            m->mem_ptr = 0;
            m->fv = NULL;
        }
    }
}
void 
q(latmat_c_copy)(latmat_c m1, latmat_c m2)
{
    assert(m1.len == m2.len);
    assert(m1.dim == m2.dim);
    assert(m1.Ls == m2.Ls);
    assert(!latmat_c_is_null(&m1));
    assert(!latmat_c_is_null(&m2));
    
    qx(vf_copy)(m1.dim, m1.Ls, m1.len, 
                m2.fv, m2.stride, m2.begin,
                m1.fv, m1.stride, m1.begin);
}
void 
q(latmat_c_swap)(latmat_c *m1, latmat_c *m2)
{
    /* laziness - the engine of progress */
    latmat_c aux;
    aux = *m1; 
    *m1 = *m2;
    *m2 = aux;
}
latmat_c
q(latmat_c_submat_col)(latmat_c m, int col, int ncol)
{
    assert(0 <= col);
    assert(0 <= ncol);
    assert(col < m.len);
    assert(col + ncol <= m.len);
    assert(!latmat_c_is_null(&(m)));
    latmat_c res;
    res.dim     = m.dim;
    res.Ls      = m.Ls;
    res.size    = m.size;
    res.begin   = m.begin + col;
    res.len     = ncol;
    res.fv      = m.fv;
    res.stride  = m.stride;
    res.mem_ptr = NULL;
    res.mem_size = 0;

    return res;
}
void
q(latmat_c_insert_col)(latmat_c m, int col, latvec_c v)
{
    assert(col < m.len);
    assert(v.dim == m.dim);
    assert(v.Ls == m.Ls);
    assert(!latmat_c_is_null(&m));
    assert(!latvec_c_is_null(&v));

    qx(vf_put)(m.dim, m.Ls,
               m.fv, m.stride, col,
               v.f);
}
void
q(latmat_c_get_col)(latmat_c m, int col, latvec_c v)
{
    assert(col < m.len);
    assert(v.dim == m.dim);
    assert(v.Ls == m.Ls);
    assert(!latmat_c_is_null(&m));
    assert(!latvec_c_is_null(&v));

    qx(vf_get)(m.dim, m.Ls,
               v.f,
               m.fv, m.stride, col);
}


/* C <- A^\dag * B, A:lat*m, B:lat*n, C:m*n */
void 
q(lat_lmH_dot_lm)(int m, int n,
               latmat_c a, 
               latmat_c b,
               doublecomplex *c, int ldc)
{
    assert(a.dim == b.dim);
    assert(a.Ls == b.Ls);
    assert(a.len == m);
    assert(b.len == n);
    assert(m <= ldc);
    assert(!latmat_c_is_null(&a));
    assert(!latmat_c_is_null(&b));
    
    int j;
    if (m == ldc) {
        memset(c, 0, m * n * sizeof(c[0]));
    } else {
        for (j = 0; j < n; j++)
            memset(c + ldc * j, 0, m * sizeof(c[0]));
    }
    
    qx(do_vfH_dot_vf)(a.dim, a.Ls,
                   (double *)c, ldc,
                   a.fv, a.stride, a.begin, a.len,
                   b.fv, b.stride, b.begin, b.len);

    if (m == ldc) {
        QMP_sum_double_array((double *)c, 2 * m * n);
    } else {
        for (j = 0; j < n; j++)
            QMP_sum_double_array((double *)(c + ldc * j), 2 * m);
    }
}
/* y <- A^\dag x, A:lat*m, x:lat, y:m */
void 
q(lat_lmH_dot_lv)(int m, 
               latmat_c a, 
               latvec_c x,  
               doublecomplex *y)
{
    assert(a.dim == x.dim);
    assert(a.Ls == x.Ls);
    assert(a.len == m);
    assert(!latmat_c_is_null(&a));
    assert(!latvec_c_is_null(&x));

    memset(y, 0, m * sizeof (y[0]));
    qx(do_vfH_dot_f)(a.dim, a.Ls,
                  (double *)y, 
                  a.fv, a.stride, a.begin, a.len,
                  x.f);
    QMP_sum_double_array((double *)y, 2 * m);
}
/* C <- A * B, A:lat*k, B:k*n, C:lat*n */
void 
q(lat_lm_dot_zm)(int n, int k, 
              latmat_c a, 
              doublecomplex *b, int ldb, 
              latmat_c c)
{
    assert(a.dim == c.dim);
    assert(a.Ls == c.Ls);
    assert(a.len == k);
    assert(c.len == n);
    assert(k <= ldb);
    assert(!latmat_c_is_null(&a));
    assert(!latmat_c_is_null(&c));

    qx(vf_dot_mz)(a.dim, a.Ls,
                  c.fv, c.stride, c.begin, c.len,
                  a.fv, a.stride, a.begin, a.len,
                  (double *)b, ldb);
}
/* y <- A * x, A:lat*n, x:n, y:lat */
void 
q(lat_lm_dot_zv)(int n, 
              latmat_c a, 
              doublecomplex *x,
              latvec_c y)
{
    assert(a.dim == y.dim);
    assert(a.Ls == y.Ls);
    assert(a.len == n);
    assert(!latmat_c_is_null(&a));
    assert(!latvec_c_is_null(&y));

    qx(vf_dot_vz)(a.dim, a.Ls,
                  y.f,
                  a.fv, a.stride, a.begin, a.len,
                  (double *)x);
}
