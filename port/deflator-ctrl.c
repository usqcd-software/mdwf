#include <assert.h>
#include <string.h>
#include <mdwf.h>

int
QX(deflator_start_load)(struct QX(Deflator) *df)
{
    if (NULL == df)
        return 1;
    if (df->loading == 1)
        return 1;
    df->loading = 1;

    return 0;
}

int
QX(deflator_stop_load)(struct QX(Deflator) *df)
{
    if (NULL == df)
        return 1;
    if (df->loading == 0)
        return 1;
    df->loading = 0;

    qx(defl_rebuild)(df);

    return 0;
}

int
QX(deflator_eigcg_resume)(struct QX(Deflator) *df)
{
    if (NULL == df)
        return 1;
    
    if ( ! df->do_eigcg)
        return q(set_error)(df->state, 0, 
                            "deflator_resume(): no EigCG workspace allocated");
    if (df->loading == 0)
        return q(set_error)(df->state, 0, 
                            "deflator_resume(): in 'loading' state");

    if (df->usize < df->umax)
        df->df_eigcg.frozen = 0;

    return 0;
}

int
QX(deflator_eigcg_stop)(struct QX(Deflator) *df)
{
    if (NULL == df)
        return 1;
    if (df->loading != 0 || ! df->do_eigcg)
        return 0;
    
    df->df_eigcg.frozen = 1;
    
    return 0;
}

int
QX(deflator_eigcg_is_stopped)(struct QX(Deflator) *df)
{
    if (df == NULL)
        return 1;
    return (!df->do_eigcg || df->df_eigcg.frozen);
}

int
QX(deflator_current_dim)(struct QX(Deflator) *df)
{
    if (df == NULL)
        return -1;
    return df->usize;
}

int
QX(deflator_eigen)(int n, double *evals,
                   struct QX(Deflator) *df)
{
    int i;

    if ((df == 0) || (evals == 0))
        return 1;

    if (df->loading)
        return q(set_error)(df->state, 0, "deflator_eigen(): in loading state");
    
    if (df->usize < n)
        n = df->usize;
    for (i = 0; i < n; i++)
        evals[i] = df->hevals[i];

    /* normal return */
    return 0; 
}

int
QX(deflator_eigcg_reset)(struct QX(Deflator) *df)
{
    if (NULL == df)
        return 1;
    if (df->loading != 0 || ! df->do_eigcg)
        return 0;

    struct qx(DeflatorEigcg) *d_e = &(df->df_eigcg);
    d_e->vsize = 0;
    memset(d_e->T, 0, d_e->vmax * d_e->vmax * sizeof(d_e->T[0]));

    return 0;
}

const char *
QX(deflator_signature)(struct QX(Deflator) *df)
{
    /* This is magic: change it if the data representation is changed. */
    return "28b0155a-9cf2-4af6-8968-5be5f80ac4d3";
}
