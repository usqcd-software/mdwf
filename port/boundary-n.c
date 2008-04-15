#include <qop-mdwf3.h>
#include "mdwf.h"

Up_project qx(up_project_n)[Q(DIM)] = {
  qx(proj_Ucg0plus),
  qx(proj_Ucg1plus),
  qx(proj_Ucg2plus),
  qx(proj_Ucg3plus)
};

Down_project qx(down_project_n)[Q(DIM)] = {
  qx(proj_g0minus),
  qx(proj_g1minus),
  qx(proj_g2minus),
  qx(proj_g3minus)
};
