#include <assert.h>
#include <mdwf.h>
#include <qmp.h>

/* replace all references to `dim' and `Ls' with these macros */
#define DEFLATOR_DIM(pstate) ((pstate)->even.full_size)
#define DEFLATOR_LS(pstate) ((pstate)->Ls)

/* allocate & free */
qx(defl_vec)
qx(defl_vec_alloc)(struct Q(State) *state)
{
    qx(defl_vec) res;
    res.dim     = DEFLATOR_DIM(state);
    res.Ls      = DEFLATOR_LS(state);
    res.mem_ptr = qx(allocate_eo)(state, &res.mem_size,
                                  (void *)&res.f, 0, 1, 0);
    if (res.mem_ptr == NULL)
        res.f = NULL;
    return res;
}
qx(defl_vec) 
qx(defl_vec_view)(struct Q(State) *state, struct Fermion *f)
{
    qx(defl_vec) res;
    res.dim     = DEFLATOR_DIM(state);
    res.Ls      = DEFLATOR_LS(state);
    res.f       = f;
    res.mem_ptr = NULL;
    res.mem_size = 0;
    return res;
}
void 
qx(defl_vec_copy)(qx(defl_vec) x, qx(defl_vec) y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!defl_vec_is_null(&x));
    assert(!defl_vec_is_null(&y));
    qx(f_copy)(y.f, y.dim, y.Ls, x.f);
}
void
qx(defl_vec_zero)(qx(defl_vec) x)
{
    assert(!defl_vec_is_null(&x));
    qx(f_zero)(x.f, x.dim, x.Ls);
}
void
qx(defl_vec_free)(struct Q(State) *state, qx(defl_vec) *v)
{
    if (NULL != v && NULL != v->mem_ptr && NULL != v->f) {
        q(free)(state, v->mem_ptr, v->mem_size);
        v->mem_ptr = NULL;
        v->f = NULL;
    }
}


int
qx(defl_vec_dotu)(
        doublecomplex *res,
        qx(defl_vec) x, 
        qx(defl_vec) y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!defl_vec_is_null(&x));
    assert(!defl_vec_is_null(&y));
    assert(NULL != res);

    double s[0];
    int flops = qx(f_dot)(s, s + 1, x.dim, x.Ls, x.f, y.f);
    QMP_sum_double_array(s, 2);
    res->r  = s[0];
    res->i  = s[1];
    
    return flops;
}
int
qx(defl_vec_scal)(double alpha, qx(defl_vec) x)
{
    assert(!defl_vec_is_null(&x));
    return qx(f_rmul1)(x.f, x.dim, x.Ls, alpha);
}

int 
qx(defl_vec_axpy)(double alpha, qx(defl_vec) x, qx(defl_vec) y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!defl_vec_is_null(&x));
    assert(!defl_vec_is_null(&y));

    return qx(f_add2)(y.f, y.dim, y.Ls, alpha, x.f);
}

int
qx(defl_vec_nrm2)(double *res, qx(defl_vec) x) 
{
    assert(!defl_vec_is_null(&x));
    assert(NULL != res);
    
    int flops = qx(f_norm)(res, x.dim, x.Ls, x.f);
    QMP_sum_double(res);
    
    return flops;
}

void /* counding is done through ws */
qx(defl_vec_linop)(
        qx(defl_vec) y, 
        qx(defl_vec) x, 
        struct MxM_workspace *ws)
{
    qx(cg_operator)(y.f, x.f, ws);
}

#if 0
/* cross-precision converters */
void 
q(latvec_zc_copy)(latvec_z x, qx(defl_vec) y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!latvec_z_is_null(&x));
    assert(!defl_vec_is_null(&y));
    /* TODO LATER copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}
void 
q(latvec_cz_copy)(qx(defl_vec) x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(x.Ls == y.Ls);
    assert(!defl_vec_is_null(&x));
    assert(!latvec_z_is_null(&y));
    /* TODO LATER copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}

#endif
    
    
qx(defl_mat) 
qx(defl_mat_alloc)(struct Q(State) *state, int ncol)
{
    qx(defl_mat) res;
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

//int
//q(latmat_convert_to_blas)(qx(defl_mat) *m) 
//{ 
//    return 0; /* stub; modify if perform any transformations */
//}
//int
//q(latmat_convert_from_blas)(qx(defl_mat) *m) 
//{ 
//    return 0; /* stub; modify if perform any transformations */
//}

qx(defl_mat) 
qx(defl_mat_view)(struct Q(State) *state, int size, struct vFermion *fv)
{
    qx(defl_mat) res;
    res.dim     = DEFLATOR_DIM(state);
    res.Ls      = DEFLATOR_LS(state);
    res.size    = size;
    res.begin   = 0;
    res.len     = size;
    res.fv      = fv;
    res.stride  = qx(strideof_vfermion)(res.dim, res.Ls);
    res.mem_ptr = NULL;
    res.mem_size= 0;
    return res;
}
void
qx(defl_mat_free)(struct Q(State) *state, qx(defl_mat) *m) 
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
qx(defl_mat_copy)(qx(defl_mat) m1, qx(defl_mat) m2)
{
    assert(m1.len == m2.len);
    assert(m1.dim == m2.dim);
    assert(m1.Ls == m2.Ls);
    assert(!defl_mat_is_null(&m1));
    assert(!defl_mat_is_null(&m2));
    
    qx(vf_copy)(m1.dim, m1.Ls, m1.len, 
                m2.fv, m2.stride, m2.begin,
                m1.fv, m1.stride, m1.begin);
}
qx(defl_mat)
qx(defl_mat_submat_col)(qx(defl_mat) m, int col, int ncol)
{
    assert(0 <= col);
    assert(0 < ncol);
    assert(col < m.len);
    assert(col + ncol <= m.len);
    assert(!defl_mat_is_null(&(m)));
    qx(defl_mat) res;
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
int
qx(defl_mat_insert_col)(qx(defl_mat) m, int col, qx(defl_vec) v)
{
    assert(col < m.len);
    assert(v.dim == m.dim);
    assert(v.Ls == m.Ls);
    assert(!defl_mat_is_null(&m));
    assert(!defl_vec_is_null(&v));

    return qx(vf_put)(m.dim, m.Ls,
                      m.fv, m.stride, col,
                      v.f);
}
int
qx(defl_mat_get_col)(qx(defl_mat) m, int col, qx(defl_vec) v)
{
    assert(col < m.len);
    assert(v.dim == m.dim);
    assert(v.Ls == m.Ls);
    assert(!defl_mat_is_null(&m));
    assert(!defl_vec_is_null(&v));

    return qx(vf_get)(m.dim, m.Ls,
                      v.f,
                      m.fv, m.stride, col);
}


/* C <- A^\dag * B, A:lat*m, B:lat*n, C:m*n */
int 
qx(defl_lmH_dot_lm)(int m, int n,
               qx(defl_mat) a, 
               qx(defl_mat) b,
               doublecomplex *c, int ldc)
{
    assert(a.dim == b.dim);
    assert(a.Ls == b.Ls);
    assert(a.len == m);
    assert(b.len == n);
    assert(m <= ldc);
    assert(!defl_mat_is_null(&a));
    assert(!defl_mat_is_null(&b));
    
    int j;
    if (m == ldc) {
        memset(c, 0, m * n * sizeof(c[0]));
    } else {
        for (j = 0; j < n; j++)
            memset(c + ldc * j, 0, m * sizeof(c[0]));
    }
    
    int flops = qx(do_vfH_dot_vf)(a.dim, a.Ls,
                                  (double *)c, ldc,
                                  a.fv, a.stride, a.begin, a.len,
                                  b.fv, b.stride, b.begin, b.len);

    if (m == ldc) {
        QMP_sum_double_array((double *)c, 2 * m * n);
    } else {
        for (j = 0; j < n; j++)
            QMP_sum_double_array((double *)(c + ldc * j), 2 * m);
    }

    return flops;
}
/* y <- A^\dag x, A:lat*m, x:lat, y:m */
int 
qx(defl_lmH_dot_lv)(int m, 
               qx(defl_mat) a, 
               qx(defl_vec) x,  
               doublecomplex *y)
{
    assert(a.dim == x.dim);
    assert(a.Ls == x.Ls);
    assert(a.len == m);
    assert(!defl_mat_is_null(&a));
    assert(!defl_vec_is_null(&x));

    memset(y, 0, m * sizeof (y[0]));
    int flops = qx(do_vfH_dot_f)(a.dim, a.Ls,
                                 (double *)y, 
                                 a.fv, a.stride, a.begin, a.len,
                                 x.f);
    QMP_sum_double_array((double *)y, 2 * m);

    return flops;
}
/* C <- A * B, A:lat*k, B:k*n, C:lat*n */
int 
qx(defl_lm_dot_zm)(int n, int k, 
              qx(defl_mat) a, 
              doublecomplex *b, int ldb, 
              qx(defl_mat) c)
{
    assert(a.dim == c.dim);
    assert(a.Ls == c.Ls);
    assert(a.len == k);
    assert(c.len == n);
    assert(k <= ldb);
    assert(!defl_mat_is_null(&a));
    assert(!defl_mat_is_null(&c));

    return qx(vf_dot_mz)(a.dim, a.Ls,
                         c.fv, c.stride, c.begin, c.len,
                         a.fv, a.stride, a.begin, a.len,
                         (double *)b, ldb);
}
/* y <- A * x, A:lat*n, x:n, y:lat */
int
qx(defl_lm_dot_zv)(int n, 
              qx(defl_mat) a, 
              doublecomplex *x,
              qx(defl_vec) y)
{
    assert(a.dim == y.dim);
    assert(a.Ls == y.Ls);
    assert(a.len == n);
    assert(!defl_mat_is_null(&a));
    assert(!defl_vec_is_null(&y));

    return qx(vf_dot_vz)(a.dim, a.Ls,
                         y.f,
                         a.fv, a.stride, a.begin, a.len,
                         (double *)x);
}
