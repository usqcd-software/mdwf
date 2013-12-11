#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>
#include <math.h>


/* XXX uses WK: df->work_c_1 (children: df->work_c_2, df->zwork) */
int
qx(defl_postamble)(
        struct QX(Deflator)         *df,
        struct qx(MxM_workspace)    *ws,
        unsigned int                options)
{
    int i;

    if (NULL == df ||
            NULL == df->state ||
            df->df_eigcg.frozen || 
            df->umax <= df->usize ||
            df->df_eigcg.vsize < df->df_eigcg.nev)
        return 0;
    struct qx(DeflatorEigcg) *d_e = &(df->df_eigcg);

    int unew = 0,
        i_v = 0;
    long int usize_old = df->usize;
    while ((df->usize < df->umax) && (i_v < d_e->nev)) {
        qx(defl_mat_get_col)(d_e->V, i_v, cur_v);
        if (qx(defl_inject_back)(df, ws, cur_v)) {
          unew ++;
        }
        i_v++;
    }
    assert(usize_old + unew == df->usize);

    if (options & QOP_MDWF_LOG_EIG_POSTAMBLE) {
        for (i = 0; i < df->usize; i++)
            qf(zprint)(df->state, "postamble", "U %4d %17.9e",
                       i, df->hevals[i]);
    }
    qx(defl_rebuild)(df);
    return unew;
}
