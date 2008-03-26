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

static unsigned int const table[] = {
    0xd6487bde, 0x9fb364c6, 0x5856f637, 0xfb376576,
    0xae83eced, 0x7e6a9a86, 0xc8d5c694, 0x097394de,
    0x210eada5, 0xdfdc1590, 0xa507d3c5, 0xabe2c7c9,
    0xa3c5d228, 0xaa1d494d, 0x2e4d8f2c, 0x68ff01f1,
    0x6d582b8d, 0x0092e7d5, 0x76ed0d80, 0xd432b139,
    0x2fda048f, 0x60932529, 0x45964081, 0x20523642,
    0x47db4bd0, 0x85b84e12, 0xdd15c8be, 0x54da7bfe,
    0x3598d90b, 0xcac5dcce, 0x7c09d671, 0xe5421d48,
    0x8587eab4, 0xaf8c30eb, 0xb9a59f46, 0xd5e916b5,
    0x096ae876, 0xe0cc2a60, 0xfe4c21e6, 0x7e6e9dc8,
    0xd07b141a, 0xb52c65d2, 0x2c1e4292, 0x4b2d5fe9,
    0x0fa492f6, 0x5cb71345, 0x0e100035, 0x9d3d1f4f,
    0xa90d97eb, 0x394aadbe, 0x9837fb8f, 0xccaeb9c3,
    0x88a3e819, 0x50fb4e2a, 0xd9cbcb32, 0xcc67f9ac,
    0x6bc4c2a7, 0x8eec0adf, 0xc44e366d, 0xc206259a,
    0xa7be9b84, 0x59a571e2, 0xc74d24f3, 0x63d4af07,
    0xc18890b2, 0x52ce5212, 0x8308ff33, 0xf2c13119,
    0xa98edfe5, 0x86caca57, 0x226c7ce7, 0xfe8fa4c1,
    0xe978cc71, 0x08a939f5, 0xf8cbff5e, 0xe23891e7,
    0x8787e33d, 0xa9aa9263, 0x8c2e426d, 0x479dc715,
    0xc15ef505, 0x6ca3b8cf, 0x91ff5b00, 0x9e9d4b3d,
    0x23ec6ea0, 0xd7b727b7, 0x63f85996, 0xfa0fbb8d,
    0xd4f0aa23, 0x2f2afc14, 0x13b49907, 0x4a495372,
    0x819a0ec8, 0xc810a2af, 0x4b4df0b6, 0x3c8d9279,
    0xd46909dd, 0xbfaf9b3e, 0xb1119c34, 0x63e34319,
    0xfe8198e9, 0x4dfedb51, 0x6f8d9723, 0x76891bf1,
    0xc6575803, 0x88062886, 0x510ed9c2, 0x08d1f204,
    0xe29ecf47, 0xf2277270, 0x4aa97600, 0x4fc3d86a,
    0x62963114, 0x1ef09c39, 0x052942ba, 0xb6de51c3,
    0xf6d2a196, 0xfab03feb, 0xd51c6ee6, 0x49769f0c,
    0xa34d55a3, 0xee198038, 0x2adad611, 0xdc895400,
    0xbb3b5779, 0x7a74bdf2, 0xbdf2a9b9, 0xb2c67331};

unsigned int
sum_init(int v)
{
    return v;
}

unsigned int
sum_add(unsigned int state, int v)
{
    return state ^ table[((state >> 5) + v + 1) % NELEM(table)];
}

double
sum_fini(unsigned int state)
{
    return state * 4.6566128730773925781e-10 - 1;
}

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

void
setup_bc(int seed)
{
    int i;
    unsigned int v;
    
    for (i = 0; i < lattice[4]; i++) {
	v = sum_init(seed);
	v = sum_add(v, i);
	v = sum_add(v, seed);
	b5[i] = sum_fini(v);
	v = sum_add(v, i);
	v = sum_add(v, seed);
	c5[i] = sum_fini(v);
    }
}

int
main(int argc, char *argv[])
{
    QMP_thread_level_t qt = QMP_THREAD_SINGLE;
    int status = 1;
    int seed_bc;
    int i;

    if (QMP_init_msg_passing(&argc, &argv, qt, &qt) != QMP_SUCCESS) {
	fprintf(stderr, "QMP_init() failed\n");
	return 1;
    }

    self = QMP_get_node_number();
    primary = QMP_is_primary_node();
    if (argc != 21) {
	zprint("21 arguments expected, found %d", argc);
	zprint("Usage: operator M_5 m seed_bc Nx Ny Nz Nt"
	       "  Lx Ly Lz Lt Ls"
	       "  x y z t s c d re/im");
	QMP_finalize_msg_passing();
	return 1;
    }
    
    M_5 = atof(argv[1]);
    m = atof(argv[2]);
    seed_bc = atoi(argv[3]);
    for (i = 0; i < 4; i++) {
	network[i] = atoi(argv[i+4]);
	lattice[i] = atoi(argv[i+8]);
	fermion_pos[i] = atoi(argv[i+13]);
    }
    lattice[4] = atoi(argv[12]);
    fermion_pos[4] = atoi(argv[17]);
    fermion_color = atoi(argv[18]);
    fermion_dirac = atoi(argv[19]);
    fermion_reim = atoi(argv[20]);

    setup_bc(seed_bc);

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
