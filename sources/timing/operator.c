#include <qmp.h>
#include <stdarg.h>
#include <qop-mdwf3.h>

#define NELEM(x) (sizeof (x) / sizeof (x[0]))

static int self;
static int primary;
static int mylattice[5];
static int mylocal[5];
static int mynetwork[4];
static int mynode[4];
static double b5[128];
static double c5[128];

xprint(char *fmt, ...)
{
  char buffer[4096];
  va_list va;

  va_start(va, fmt);
  vsprintf(buffer, fmt, va);
  va_end(va);
  printf("[%04d] %s\n", self, buffer);
  fflush(stdout);
}

static void
zprint(char *fmt, ...)
{
  va_list va;
  if (primary == 0)
    return;

  va_start(va, fmt);
  vprintf(fmt, va);
  va_end(va);
}

static void
xshowv(char *name, const int v[4])
{
  xprint("%s: %d %d %d %d", name, v[0], v[1], v[2], v[3]);
}

static void
zshowv4(char *name, const int v[4])
{
  zprint("%s: %d %d %d %d\n", name, v[0], v[1], v[2], v[3]);
}

static void
zshowv5(char *name, const int v[5])
{
  zprint("%s: %d %d %d %d %d\n", name, v[0], v[1], v[2], v[3], v[4]);
}

static void
getsub(int lo[4], int hi[4], const int node[4], void *env)
{
  int i;
  for (i = 0; i < 4; i++) {
    lo[i] = (mylattice[i] * node[i]) / mynetwork[i];
    hi[i] = (mylattice[i] * (node[i] + 1)) / mynetwork[i];
  }
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

static double
read_gauge(int dir,
	   const int pos[4],
	   int a,
	   int b,
	   int re_im,
	   void *env)
{
  return 1.0;
}

static double
read_fermion(const int pos[5],
	     int a,
	     int b,
	     int re_im,
	     void *env)
{
  return 1.0;
}

static void
report_perf(const char *name, struct QOP_MDWF_State *state)
{
  double sec;
  long long flops;
  long long bytes;
  
  QOP_MDWF_performance(&sec, &flops, &bytes, state);
  zprint("perf(%s): %g MFlops/sec, %g MBytes/sec [%g %lld %lld]\n",
	 name, 1e-6 * flops/sec, 1e-6 * bytes/sec,
	 sec, flops, bytes);
}

static int
do_f3(struct QOP_MDWF_State *state, struct QOP_MDWF_Parameters *params)
{
  struct QOP_F3_MDWF_Gauge *U;
  struct QOP_F3_MDWF_Fermion *F;
  struct QOP_F3_MDWF_Fermion *R;

  if (QOP_F3_MDWF_import_gauge(&U, state, read_gauge, NULL)) {
    zprint("f3: import gauge failed\n");
    goto no_U;
  }
  if (QOP_F3_MDWF_import_fermion(&F, state, read_fermion, NULL)) {
    zprint("f3: import fermion failed\n");
    goto no_F;
  }
  if (QOP_F3_MDWF_allocate_fermion(&R, state)) {
    zprint("f3: allocate fermion failed\n");
    goto no_R;
  }

  if (QOP_F3_MDWF_DDW_operator(R, params, U, F)) {
    zprint("f3: operator failed: %s", QOP_MDWF_error(state));
    goto no_X;
  }
  report_perf("f3", state);

  QOP_F3_MDWF_free_fermion(&R);    
  QOP_F3_MDWF_free_fermion(&F);    
  QOP_F3_MDWF_free_gauge(&U);
  return 0;
 no_X:
  QOP_F3_MDWF_free_fermion(&R);    
 no_R:
  QOP_F3_MDWF_free_fermion(&F);    
 no_F:
  QOP_F3_MDWF_free_gauge(&U);  
 no_U:
  return 1;
}

static int
do_d3(struct QOP_MDWF_State *state, struct QOP_MDWF_Parameters *params)
{
  struct QOP_D3_MDWF_Gauge *U;
  struct QOP_D3_MDWF_Fermion *F;
  struct QOP_D3_MDWF_Fermion *R;

  if (QOP_D3_MDWF_import_gauge(&U, state, read_gauge, NULL)) {
    zprint("d3: import gauge failed\n");
    goto no_U;
  }
  if (QOP_D3_MDWF_import_fermion(&F, state, read_fermion, NULL)) {
    zprint("d3: import fermion failed\n");
    goto no_F;
  }
  if (QOP_D3_MDWF_allocate_fermion(&R, state)) {
    zprint("d3: allocate fermion failed\n");
    goto no_R;
  }

  QOP_D3_MDWF_DDW_operator(R, params, U, F);
  report_perf("d3", state);

  QOP_D3_MDWF_free_fermion(&R);    
  QOP_D3_MDWF_free_fermion(&F);    
  QOP_D3_MDWF_free_gauge(&U);
  return 0;
 no_R:
  QOP_D3_MDWF_free_fermion(&F);    
 no_F:
  QOP_D3_MDWF_free_gauge(&U);  
 no_U:
  return 1;
}

int
main(int argc, char *argv[])
{
  struct QOP_MDWF_State *mdwf_state = NULL;
  struct QOP_MDWF_Parameters *mdwf_params = NULL;
  QMP_thread_level_t qt = QMP_THREAD_SINGLE;
  int status = 1;
  int i;

  if (QMP_init_msg_passing(&argc, &argv, qt, &qt) != QMP_SUCCESS) {
    fprintf(stderr, "QMP_init() failed\n");
    return 1;
  }

  for (i = 0; i < NELEM(b5); i++) {
    b5[i] = 0.1 * i * (NELEM(b5) - i);
    c5[i] = 0.1 * i * i * (NELEM(b5) - i);
  }

  self = QMP_get_node_number();
  primary = QMP_is_primary_node();
  for (i = 0; i < argc; i++)
    zprint("arg[%d]=%s\n", i, argv[i]);
  if (argc != 10) {
    zprint("10 arguments expected, found %d\n", argc);
    QMP_finalize_msg_passing();
    return 1;
  }

  for (i = 0; i < 4; i++) {
    mynetwork[i] = atoi(argv[i+1]);
    mylocal[i] = atoi(argv[i+5]);
    mylattice[i] = mylocal[i] * mynetwork[i];
  }
  mylocal[4] = mylattice[4] = atoi(argv[9]);

  zshowv4("network", mynetwork);
  zshowv5("local lattice", mylocal);
  zshowv5("lattice", mylattice);

  if (QMP_declare_logical_topology(mynetwork, 4) != QMP_SUCCESS) {
    zprint("declare_logical_top failed\n");
    goto end;
  }

  getv(mynode, 0, QMP_get_logical_number_of_dimensions(),
       QMP_get_logical_coordinates());

  xshowv("node", mynode);

  if (QOP_MDWF_init(&mdwf_state,
		    mylattice, mynetwork, mynode, primary,
		    getsub, NULL)) {
    zprint("MDWF_init() failed\n");
    goto end;
  }

  zprint("MDWF_init() done\n");

  if (QOP_MDWF_set_generic(&mdwf_params, mdwf_state, b5, c5, 0.123, 0.05)) {
    zprint("MDW_set_generic() failed\n");
    goto end;
  }
  zprint("MDWF_set_generic() done\n");

  if (do_f3(mdwf_state, mdwf_params)) {
    zprint("float test failed\n");
    goto end;
  }
  if (do_d3(mdwf_state, mdwf_params)) {
    zprint("double test failed\n");
    goto end;
  }

  QOP_MDWF_fini(&mdwf_state);

  status = 0;
 end:
  QMP_finalize_msg_passing();
  return status;
}
