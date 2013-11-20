#include <mdwf.h>

int
QX(MxM_poly)(struct QX(HalfFermion) *result,
             struct QX(HalfFermion) *result_prev,
             const struct Q(Parameters) *params,
             const struct QX(Gauge) *gauge,
             const struct QX(HalfFermion) *fermion,
             int poly_n, 
             const double poly_a[], 
             const double poly_b[], 
             const double poly_c[])
{
    DECLARE_STATE;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    size_t alloc_size;
    void *alloc_ptr;
    void *aligned_ptr;
    struct Fermion *tmp_o;
    struct Fermion *tmp_e1,
                   *tmp_e2;
    int Ls;
    int e_size;
    int i;
    struct Fermion *poly_ws_e[3];
    struct Fermion *poly_v0_e, 
                   *poly_v1_e, 
                   *poly_v2_e,
                   *aux;
    
    CHECK_ARG0(result);
    CHECK_ARGn(params, "MxM_poly");
    CHECK_ARGn(gauge, "MxM_poly");
    CHECK_ARGn(fermion, "MxM_poly");

    Ls      = state->Ls;
    e_size  = state->even.full_size;

    if (poly_n <= 0) 
        return q(set_error)(state, 0,
                            "MxM_poly(): degree must be > 0");

    if (0 == poly_a 
            || 0 == poly_b
            || 0 == poly_c)
        return q(set_error)(state, 0,
                            "MxM_poly(): NULL coefficient tables");
    
    if (q(setup_comm)(state, sizeof (REAL)))
        return q(set_error)(state, 0,
                            "MxM_poly(): communication setup failed");

    /* ws for 3-term recurrence & 1 even, 1 odd tmp vectors M^\dag.M */
    alloc_ptr = qx(allocate_eo)(state, &alloc_size, &aligned_ptr, 0, 4, 1);
    if (alloc_ptr == 0)
        return q(set_error)(state, 0, "MxM_poly(): not enough memory");

    poly_ws_e[0]    = result->even;
    poly_ws_e[1]    = aligned_ptr;
    poly_ws_e[2]    = aligned_ptr = qx(step_even)(state, aligned_ptr);
    tmp_e1          = aligned_ptr = qx(step_even)(state, aligned_ptr);
    tmp_e2          = aligned_ptr = qx(step_even)(state, aligned_ptr);
    tmp_o           = aligned_ptr = qx(step_even)(state, aligned_ptr);
    /* arrange pointers so that (poly_n-1) cyc.shifts v0<-v1<-v2
       set v1 == result->even
       originally, orig{v}[k] <- ws[(-k+n)%3]
       after (n-1) shifts, v[k] <- orig{v}[(k+n-1)%3] = ws[(1-k)%3], 
       so that v[1] <- ws[0] = result->even */
    poly_v0_e       = poly_ws_e[(0 + poly_n) % 3];
    poly_v1_e       = poly_ws_e[(2 + poly_n) % 3];
    poly_v2_e       = poly_ws_e[(1 + poly_n) % 3];
   
    BEGIN_TIMING(state); 
    /* v0 <- c0 * v */
    qx(f_copy) (poly_v0_e, e_size, Ls, fermion->even); 
    qx(f_rmul1)(poly_v0_e, e_size, Ls, poly_c[0]);
    /* v1 <- (a0 + b0*A) . v */
    qx(op_even_M) (tmp_e1, state, params, gauge->data, fermion->even,
                  &flops, &sent, &received, tmp_o);
    qx(op_even_Mx)(poly_v1_e, state, params, gauge->data, tmp_e1,
                  &flops, &sent, &received, tmp_e2, tmp_o);
    qx(f_rmul1)(poly_v1_e, e_size, Ls, poly_b[0]);
    qx(f_add2) (poly_v1_e, e_size, Ls, poly_a[0], fermion->even);
    
    for (i = 1; i < poly_n ; i++) {
        /* v2 = A.v1 */
        qx(op_even_M)(tmp_e1, state, params, gauge->data, poly_v1_e,
                      &flops, &sent, &received, tmp_o);
        qx(op_even_Mx)(poly_v2_e, state, params, gauge->data, tmp_e1,
                      &flops, &sent, &received, tmp_e2, tmp_o);
        /* v2 *= b */
        qx(f_rmul1)(poly_v2_e, e_size, Ls, poly_b[i]);
        /* v2 += a*v1 */
        qx(f_add2) (poly_v2_e, e_size, Ls, poly_a[i], poly_v1_e);
        /* v2 += c*v0 */
        qx(f_add2) (poly_v2_e, e_size, Ls, poly_c[i], poly_v0_e);

        /* shift v0 <- v1 <- v2*/
        aux         = poly_v0_e;
        poly_v0_e   = poly_v1_e;
        poly_v1_e   = poly_v2_e;
        poly_v2_e   = aux;
        /* the last computed poly is in v1 */
    }
    END_TIMING(state, flops, sent, received);

    if (NULL != result_prev) {
        CHECK_ARGn(result_prev, "MxM_poly");
        qx(f_copy)(result_prev->even, e_size, Ls, poly_v0_e);
    }

    q(free)(state, alloc_ptr, alloc_size);

    return 0;
} 

