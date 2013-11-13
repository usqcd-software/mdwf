#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

#if QOP_MDWF_DEFAULT_PRECISION == 'F'
int
QX(deflator_extract_vector)(struct QX(HalfFermion) *hf,
                            const struct Q(Deflator) *d,
                            int idx)
{
  DECLARE_STATE;
  latvec_c vf;

  CHECK_ARG0(hf);
  CHECK_ARGn(d, "deflator_extract_vector");

  if (idx < 0 || idx >= d->usize)
    return q(set_error)(state, 0, "deflator_extract_vector(): index out of range");

  vf = q(latvec_c_view)(d->dim, d->Ls, hf->even);
  q(latmat_c_get_col)(d->U, idx, vf);

  return 0;
}
#endif /* QOP_MDWF_DEFAULT_PRECISION == 'F' */
