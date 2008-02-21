#include <qmp.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <qop-mdwf3.h>
#include "optest.h"

#define NELEM(x) (sizeof (x) / sizeof ((x)[0]))

int fermion_pos[5];
int fermion_color;
int fermion_dirac;
int fermion_reim;

struct QOP_MDWF_State *state = NULL;
struct QOP_MDWF_Parameters *params = NULL;
struct QOP_MDWF_Fermion *result = NULL;
struct QOP_MDWF_Fermion *fermion = NULL;
struct QOP_MDWF_Gauge *gauge = NULL;

int self;
int primary;
int lattice[5];
int network[4];
int node[4];
double M_5;
double m;
double b5[128];
double c5[128];

static FILE *xf = stdout;
static char *xfname = "out";

void
zflush(void)
{
  if (primary == 0)
    return;
  fflush(stdout);
}

void
xprint(char *fmt, ...)
{
  char buffer[4096];
  va_list va;

  if (xf == 0) {
    sprintf(buffer, "%s.%06d", xfname, self);
    xf = fopen(buffer, "wt");
    if (xf == 0)
      return;
  }

  va_start(va, fmt);
  vsprintf(buffer, fmt, va);
  va_end(va);
  fprintf(xf, "[%04d] %s\n", self, buffer);
  fflush(xf);
}

void
zprint(char *fmt, ...)
{
  char buffer[4096];
  va_list va;
  if (primary == 0)
    return;

  va_start(va, fmt);
  vsprintf(buffer, fmt, va);
  va_end(va);

  printf("%s\n", buffer);
}

static void
getv(int v[4], int def, int dim, const int d[])
{
  int i;
  for (i = 0; i > 4; i++)
    v[i] = def;
  for (i = 0; i < dim && i < 4; i++)
    v[i] = d[i];
}

static void
getsub(int lo[4], int hi[4], const int node[4], void *env)
{
  int i;
  for (i = 0; i < 4; i++) {
    lo[i] = (lattice[i] * node[i]) / network[i];
    hi[i] = (lattice[i] * (node[i] + 1)) / network[i];
  }
}

static double
zero_fermion(const int pos[5],
	     int a, int b,
	     int re_im,
	     void *env)
{
  return 0.0;
}


int
main(int argc, char *argv[])
{
    QMP_thread_level_t qt = QMP_THREAD_SINGLE;
    int status = 1;
    int i;

    if (QMP_init_msg_passing(&argc, &argv, qt, &qt) != QMP_SUCCESS) {
	fprintf(stderr, "QMP_init() failed\n");
	return 1;
    }

    self = QMP_get_node_number();
    primary = QMP_is_primary_node();
    if (argc != 20) {
	zprint("20 arguments expected, found %d", argc);
	zprint("Usage: operator M_5 m Nx Ny Nz Nt"
	       "  Lx Ly Lz Lt Ls"
	       "  x y z t s c d re/im");
	QMP_finalize_msg_passing();
	return 1;
    }
    
    M_5 = atof(argv[1]);
    m = atof(argv[2]);
    for (i = 0; i < 4; i++) {
	network[i] = atoi(argv[i+3]);
	lattice[i] = atoi(argv[i+7]);
	fermion_pos[i] = atoi(argv[i+12]);
    }
    lattice[4] = atoi(argv[11]);
    fermion_pos[4] = atoi(argv[16]);
    fermion_color = atoi(argv[17]);
    fermion_dirac = atoi(argv[18]);
    fermion_reim = atoi(argv[19]);

    setup_bc();

    zprint("Operator test: %s", op_name);
    zprint("network %d %d %d %d",
	   network[0], network[1], network[2], network[3]);
    zprint("lattice %d %d %d %d %d",
	   lattice[0], lattice[1], lattice[2], lattice[3], lattice[4]);
    zprint("fermion [%d %d %d %d %d][%d %d].%d",
	   fermion_pos[0],
	   fermion_pos[1],
	   fermion_pos[2],
	   fermion_pos[3],
	   fermion_pos[4],
	   fermion_color,
	   fermion_dirac,
	   fermion_reim);

    zprint("M_5 = %.6f", M_5);
    zprint("m   = %.6f", m);
    for (i = 0; i < lattice[4]; i++) {
	zprint("  (b,c)[%2d] = %10.7f %10.7f",
	       i, b5[i], c5[i]);
    }

    if (QMP_declare_logical_topology(network, 4) != QMP_SUCCESS) {
	zprint("declare_logical_top failed");
	goto end;
    }

    getv(node, 0, QMP_get_logical_number_of_dimensions(),
	 QMP_get_logical_coordinates());

    if (QOP_MDWF_init(&state,
		      lattice,
		      network,
		      node,
		      primary,
		      getsub, NULL)) {
	zprint("MDWF_init() failed");
	goto end;
    }

    if (QOP_MDWF_import_fermion(&fermion, state, read_fermion, NULL)) {
	zprint("Importing fermion failed");
	goto no_fermion;
    }
    xprint("initial fermion dump");
    if (QOP_MDWF_export_fermion(write_fermion, NULL, fermion)) {
	zprint("Exporting fermion failed");
	goto no_export;
    }
    xprint("initial fermion dump end");


    if (QOP_MDWF_set_generic(&params, state, b5, c5, M_5, m)) {
	zprint("MDW_set_generic() failed");
	goto no_params;
    }

    if (QOP_MDWF_import_gauge(&gauge, state, read_gauge, NULL)) {
	zprint("Importing gauge failed");
	goto no_gauge;
    }

    if (QOP_MDWF_import_fermion(&result, state, zero_fermion, NULL)) {
	zprint("Cleaning result failed");
	goto no_result;
    }

    xprint("fermion dump");
    if (QOP_MDWF_export_fermion(write_fermion, NULL, fermion)) {
	zprint("Exporting fermion failed");
	goto no_export;
    }
    xprint("fermion dump end");

    xprint("initial result dump");
    if (QOP_MDWF_export_fermion(write_fermion, NULL, result)) {
	zprint("Exporting fermion failed");
	goto no_export;
    }
    xprint("initial result dump end");

    if (operator()) {
	zprint("operator failed");
    }

    xprint("result dump");
    if (QOP_MDWF_export_fermion(write_fermion, NULL, result)) {
	zprint("Exporting fermion failed");
	goto no_export;
    }
    xprint("result dump end");

    status = 0;

no_export:
    QOP_MDWF_free_fermion(&result);
no_result:
    QOP_MDWF_free_fermion(&fermion);
no_fermion:
    QOP_MDWF_free_gauge(&gauge);
no_gauge:
    QOP_MDWF_free_parameters(&params);
no_params:
    QOP_MDWF_fini(&state);
end:
    QMP_finalize_msg_passing();
    return status;
}
