#include <qop-mdwf3.h>
#include "../../port/mdwf.h"
#include "op-routines.h"

Up_project up_project_x[Q(DIM)] = {
  qx(proj_Ucg0minus),
  qx(proj_Ucg1minus),
  qx(proj_Ucg2minus),
  qx(proj_Ucg3minus)
};

Down_project down_project_x[Q(DIM)] = {
  qx(proj_g0plus),
  qx(proj_g1plus),
  qx(proj_g2plus),
  qx(proj_g3plus)
};

