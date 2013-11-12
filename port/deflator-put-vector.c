#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <mdwf.h>

#if QOP_MDWF_DEFAULT_PRECISION == 'F'
int
QX(deflator_add_vector)(struct Q(Deflator) *delf_ptr,
                        const struct QX(HalfFermion) *hfermion_ptr)
{
  /* XXX */
  return 1;
}
#endif /* QOP_MDWF_DEFAULT_PRECISION == 'F' */
