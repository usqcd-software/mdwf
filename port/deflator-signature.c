#include <mdwf.h>

const char *
Q(deflator_signature)(struct Q(Deflator) *d)
{
  /* This is magic: change it if the data representation is changed. */
  return "28b0155a-9cf2-4af6-8968-5be5f80ac4d3";
}
