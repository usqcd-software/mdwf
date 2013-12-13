#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <math.h>
#include <assert.h>
#include <mdwf.h>


/* XXX uses WK: df->work_c_1 (children: df->work_c_2, df->zwork) */
int
qx(defl_eigcg_postamble)(
        struct QX(Deflator)         *df,
        struct qx(MxM_workspace)    *ws,
        unsigned int                options)
{
    int i;
    long long fl = 0;

    if (NULL == df || NULL == df->state) 
        return 1;
    if (df->df_eigcg.frozen
            || df->umax <= df->usize
            || df->df_eigcg.vsize < df->df_eigcg.nev)
        return 0;
    struct qx(DeflatorEigcg) *d_e = &(df->df_eigcg);

    int unew = 0,
        i_v = 0;
    long int usize_old = df->usize;
    qx(defl_vec) cur_v = df->work_c_1;
    while ((df->usize < df->umax) && (i_v < d_e->nev)) {
        fl += qx(defl_mat_get_col)(d_e->V, i_v, cur_v);
        if (qx(defl_inject_back)(df, ws, cur_v)) {
          unew ++;
        }
        i_v++;
    }
    assert(usize_old + unew == df->usize);

    if (options & QOP_MDWF_LOG_EIG_POSTAMBLE) {
        for (i = 0; i < df->usize; i++)
            qf(zprint)(df->state, "defl_eigcg_postamble", "U %4d %17.9e",
                       i, df->hevals[i]);
    }
    qx(defl_rebuild)(df);

    df->state->flops += fl;
    return unew;
}
