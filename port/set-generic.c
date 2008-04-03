#include <mdwf.h>

int
Q(set_generic)(struct Q(Parameters) **params_ptr,
	       struct Q(State) *state,
	       const double b_5[],
	       const double c_5[],
	       double M_5,
	       double m)
{
  char *msg = NULL;
  int i;
  int Ls;
  struct Q(Parameters) *params;
  int abs;
  int iabs;

#define CHECK(p, m) do { if ((p) == NULL) { msg = (m); goto error; }} while(0)
#define u_a(j) (c_5[j] * (M_5 + 4) - 1)
#define v_a(j) (-m * u_a(j))
#define w_a(j) (b_5[j] * (M_5 + 4) + 1)
#define u_b(j) (c_5[j])
#define v_b(j) (-m * u_b(j))
#define w_b(j) (b_5[j])

  if (state == NULL || state->error_latched)
    return 1;

  Ls = state->Ls;
  if (params_ptr == NULL)
    return q(set_error)(state, 0, "set_generic(): NULL pointer");

  *params_ptr = NULL;
  params = q(malloc)(state, sizeof (struct Q(Parameters)));
  CHECK(params, "set_generic(): Not enough memory for parameters");
  memset(params, 0, sizeof (struct Q(Parameters)));

  /* size of A & B tables */
  abs = q(sizeof_ABTable)(Ls);
  /* A tables */
  params->ApTable = q(malloc)(state, abs);
  CHECK(params->ApTable, "set_generic(): Not enough memory for A/+ table");
  q(put_ABTable)(params->ApTable, 0, w_a(0), v_a(0));
  for (i = 1; i < Ls; i++)
      q(put_ABTable)(params->ApTable, i, w_a(i), u_a(i));
  params->AmTable = q(malloc)(state, abs);
  CHECK(params->AmTable, "set_generic(): Not enough memory for A/- table");
  q(put_ABTable)(params->AmTable, 0, w_a(Ls-1), v_a(Ls-1));
  for (i = 1; i < Ls; i++)
      q(put_ABTable)(params->AmTable, i, w_a(Ls-1-i), u_a(Ls-1-i));

  /* B tables */
  params->BpTable = q(malloc)(state, abs);
  CHECK(params->BpTable, "set_generic(): Not enough memory for B/+ table");
  q(put_ABTable)(params->BpTable, 0, w_b(0), v_b(0));
  for (i = 1; i < Ls; i++)
      q(put_ABTable)(params->BpTable, i, w_b(i), u_b(i));
  params->BmTable = q(malloc)(state, abs);
  CHECK(params->BmTable, "set_generic(): Not enough memory for B/- table");
  q(put_ABTable)(params->BmTable, 0, w_b(Ls-1), v_b(Ls-1));
  for (i = 1; i < Ls; i++)
      q(put_ABTable)(params->BmTable, i, w_b(Ls-1-i), u_b(Ls-1-i));

  /* A* tables */
  params->AxpTable = q(malloc)(state, abs);
  CHECK(params->AxpTable, "set_generic(): Not enough memory for A*/+ table");
  q(put_ABTable)(params->AxpTable, 0, w_a(0), v_a(Ls-1));
  for (i = 1; i < Ls; i++)
      q(put_ABTable)(params->AxpTable, i, w_a(i), u_a(i-1));
  params->AxmTable = q(malloc)(state, abs);
  CHECK(params->AxmTable, "set_generic(): Not enough memory for A*/- table");
  q(put_ABTable)(params->AxmTable, 0, w_a(Ls-1), v_a(0));
  for (i = 1; i < Ls; i++)
      q(put_ABTable)(params->AxmTable, i, w_a(Ls-1-i), u_a(Ls-i));

  /* B* tables */
  params->BxpTable = q(malloc)(state, abs);
  CHECK(params->BxpTable, "set_generic(): Not enough memory for B*/+ table");
  q(put_ABTable)(params->BxpTable, 0, w_b(0), v_b(Ls-1));
  for (i = 1; i < Ls; i++)
      q(put_ABTable)(params->BxpTable, i, w_b(i), u_b(i-1));
  params->BxmTable = q(malloc)(state, abs);
  CHECK(params->BxmTable, "set_generic(): Not enough memory for B*/- table");
  q(put_ABTable)(params->BxmTable, 0, w_b(Ls-1), v_b(0));
  for (i = 1; i < Ls; i++)
      q(put_ABTable)(params->BxmTable, i, w_b(Ls-1-i), u_b(Ls-i));

  /* sizeof 1/A tables */
  iabs = q(sizeof_ABiTable)(Ls);
  /* 1/A tables */
  params->AipTable = q(malloc)(state, iabs);
  CHECK(params->AipTable, "set_generic(): Not enough memory for 1/A + table");
  {
      double ak = - v_a(0) / w_a(Ls-1);
      for (i = Ls; --i;) {
	  double bk = - u_a(i) / w_a(i);
	  double ck = 1.0 / w_a(i);
	  q(put_ABiTable)(params->AipTable, i, ak, bk, ck);
	  ak = - ak * u_a(i) / w_a(i-1);
      }
      q(put_ABiTableZ)(params->AipTable, 1.0 / (w_a(0) * (1.0 - ak)));
  }
  params->AimTable = q(malloc)(state, iabs);
  CHECK(params->AimTable, "set_generic(): Not enough memory for 1/A - table");
  {
      double ak = - v_a(Ls-1) / w_a(0);
      for (i = 0; i < Ls - 1; i++) {
	  double bk = -u_a(i) / w_a(i);
	  double ck = 1.0 / w_a(i);
	  q(put_ABiTable)(params->AimTable, i + 1, ak, bk, ck);
	  ak = - ak * u_a(i) / w_a(i+1);
      }
      q(put_ABiTableZ)(params->AimTable, 1.0 / (w_a(Ls-1) * (1.0 - ak)));
  }

  /* 1/B tables */
  params->BipTable = q(malloc)(state, iabs);
  CHECK(params->BipTable, "set_generic(): Not enough memory for 1/B + table");
  {
      double ak = - v_b(0) / w_b(Ls-1);
      for (i = Ls; --i;) {
	  double bk = - u_b(i) / w_b(i);
	  double ck = 1.0 / w_b(i);
	  q(put_ABiTable)(params->BipTable, i, ak, bk, ck);
	  ak = - ak * u_b(i) / w_b(i-1);
      }
      q(put_ABiTableZ)(params->BipTable, 1.0 / (w_b(0) * (1.0 - ak)));
  }
  params->BimTable = q(malloc)(state, iabs);
  CHECK(params->BimTable, "set_generic(): Not enough memory for 1/B - table");
  {
      double ak = - v_b(Ls-1) / w_b(0);
      for (i = 0; i < Ls - 1; i++) {
	  double bk = -u_b(i) / w_b(i);
	  double ck = 1.0 / w_b(i);
	  q(put_ABiTable)(params->BimTable, i + 1, ak, bk, ck);
	  ak = - ak * u_b(i) / w_b(i+1);
      }
      q(put_ABiTableZ)(params->BimTable, 1.0 / (w_b(Ls-1) * (1.0 - ak)));
  }

  /* 1/A+ tables */
  params->AxipTable = q(malloc)(state, iabs);
  CHECK(params->AxipTable, "set_generic(): Not enough memory for 1/A* + table");
  {
      double ak = - v_a(0) / w_a(0);
      for (i = 1; i < Ls; i++) {
	  double bk = -u_a(i) / w_a(i-1);
	  double ck = 1.0 / w_a(i-1);
	  q(put_ABiTable)(params->AxipTable, i, ak, bk, ck);
	  ak = - ak * u_a(i) / w_a(i);
      }
      q(put_ABiTableZ)(params->AxipTable, 1.0 / (w_a(Ls-1) * (1.0 - ak)));
  }
  params->AximTable = q(malloc)(state, iabs);
  CHECK(params->AximTable, "set_generic(): Not enough memory for 1/A* - table");
  {
      double ak = - v_a(Ls-1) / w_a(Ls-1);
      for (i = Ls; --i;) {
	  double bk = - u_a(i-1) / w_a(i);
	  double ck = 1.0 / w_a(i);
	  q(put_ABiTable)(params->AximTable, i, ak, bk, ck);
	  ak = - ak * u_a(i-1) / w_a(i-1);
      }
      q(put_ABiTableZ)(params->AximTable, 1.0 / (w_a(0) * (1.0 - ak)));
  }

  /* 1/B+ tables */
  params->BxipTable = q(malloc)(state, iabs);
  CHECK(params->BxipTable, "set_generic(): Not enough memory for 1/B* + table");
  {
      double ak = - v_b(0) / w_b(0);
      for (i = 1; i < Ls; i++) {
	  double bk = -u_b(i) / w_b(i-1);
	  double ck = 1.0 / w_b(i-1);
	  q(put_ABiTable)(params->BxipTable, i, ak, bk, ck);
	  ak = - ak * u_b(i) / w_b(i);
      }
      q(put_ABiTableZ)(params->BxipTable, 1.0 / (w_b(Ls-1) * (1.0 - ak)));
  }
  params->BximTable = q(malloc)(state, iabs);
  CHECK(params->BximTable, "set_generic(): Not enough memory for 1/B* - table");
  {
      double ak = - v_b(Ls-1) / w_b(Ls-1);
      for (i = Ls; --i;) {
	  double bk = - u_b(i-1) / w_b(i);
	  double ck = 1.0 / w_b(i);
	  q(put_ABiTable)(params->BximTable, i, ak, bk, ck);
	  ak = - ak * u_b(i-1) / w_b(i-1);
      }
      q(put_ABiTableZ)(params->BximTable, 1.0 / (w_b(0) * (1.0 - ak)));
  }

  BEGIN_TIMING(state);
  params->state = state;
  *params_ptr = params;
  END_TIMING(state, 0, 0, 0);

  return 0;
#undef u_a
#undef v_a
#undef w_a
#undef u_b
#undef v_b
#undef w_b
#undef CHECK
 error:
  Q(free_parameters)(&params);
  return q(set_error)(state, 0, msg);
}
