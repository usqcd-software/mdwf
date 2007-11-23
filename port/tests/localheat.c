#include <qmp.h>
#include <stdlib.h>
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

#define SHOW_INTERVAL 60.0
static double total_sec;

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
  FILE *xf = stderr;
  va_list va;

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

void
xshowv(char *name, const int v[4])
{
  xprint("%s: %d %d %d %d", name, v[0], v[1], v[2], v[3]);
}

static void
zshowv4(char *name, const int v[4])
{
  zprint("%s: %d %d %d %d", name, v[0], v[1], v[2], v[3]);
}

static void
zshowv5(char *name, const int v[5])
{
  zprint("%s: %d %d %d %d %d", name, v[0], v[1], v[2], v[3], v[4]);
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

static double t_sec;
static long long t_flops;
static long long t_send;
static long long t_receive;
static int rounds;

static void
start_perf(void)
{
  t_sec = 0.0;
  t_flops = 0;
  t_send = 0;
  t_receive = 0;
  rounds = 0;
}
  
static int
report_perf(const char prec,
	    struct QOP_MDWF_State *state,
	    double *run_time,
	    double *report_moment)
{
  double sec;
  long long flops;
  long long send;
  long long receive;
  
  QOP_MDWF_performance(&sec, &flops, &send, &receive, state);
  t_sec += sec;
  t_flops += flops;
  t_send += send;
  t_receive += receive;
  rounds++;

  *run_time += sec;
  if (*run_time < *report_moment)
    return 1;

  *report_moment += SHOW_INTERVAL;

  zprint("perf(%c3 @ %d): %g MFlops/sec, snd %g MBytes/sec, rcv %g MBytes/sec "
	 "[%g %lld %lld %lld]",
	 prec, rounds,
	 1e-6 * t_flops/t_sec,
	 1e-6 * t_send/t_sec,
	 1e-6 * t_receive/t_sec,
	 t_sec, t_flops, t_send, t_receive);
  zflush();
  if (*run_time < total_sec)
    return 1;

  return 0;
}

static int
do_run(struct QOP_MDWF_State *state, struct QOP_MDWF_Parameters *params)
{
  double run_time = 0.0;
  double next_report = SHOW_INTERVAL;
  struct QOP_MDWF_Gauge *U;
  struct QOP_MDWF_Fermion *F;
  struct QOP_MDWF_Fermion *R;

  if (QOP_MDWF_import_gauge(&U, state, read_gauge, NULL)) {
    zprint("%c3: import gauge failed", QOP_MDWF_DEFAULT_PRECISION);
    goto no_U;
  }
  if (QOP_MDWF_import_fermion(&F, state, read_fermion, NULL)) {
    zprint("%c3: import fermion failed", QOP_MDWF_DEFAULT_PRECISION);
    goto no_F;
  }
  if (QOP_MDWF_allocate_fermion(&R, state)) {
    zprint("%c3: allocate fermion failed", QOP_MDWF_DEFAULT_PRECISION);
    goto no_R;
  }

  start_perf();
  do {
    if (QOP_MDWF_DDW_operator(R, params, U, F)) {
      xprint("%c3: operator failed: %s", QOP_MDWF_error(state),
	     QOP_MDWF_DEFAULT_PRECISION);
      goto no_X;
    }
  } while(report_perf(QOP_MDWF_DEFAULT_PRECISION, state,
		      &run_time, &next_report));

  QOP_MDWF_free_fermion(&R);    
  QOP_MDWF_free_fermion(&F);    
  QOP_MDWF_free_gauge(&U);
  return 0;
 no_X:
  QOP_MDWF_free_fermion(&R);    
 no_R:
  QOP_MDWF_free_fermion(&F);    
 no_F:
  QOP_MDWF_free_gauge(&U);  
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
  if (argc != 7) {
    zprint("7 arguments expected, found %d", argc);
    zprint("usage: localheat Lx Ly Lz Lt Ls time");
    QMP_finalize_msg_passing();
    return 1;
  }

  for (i = 0; i < 4; i++) {
    mynetwork[i] = 1;
    mylocal[i] = atoi(argv[i+1]);
    mylattice[i] = mylocal[i] * mynetwork[i];
  }
  mylocal[4] = mylattice[4] = atoi(argv[5]);
  total_sec = atoi(argv[6]);

  zshowv4("network", mynetwork);
  zshowv5("local lattice", mylocal);
  zshowv5("lattice", mylattice);
  zprint("total requested runtime %.0f sec", total_sec);

#if 0
  if (QMP_declare_logical_topology(mynetwork, 4) != QMP_SUCCESS) {
    zprint("declare_logical_top failed");
    goto end;
  }

  getv(mynode, 0, QMP_get_logical_number_of_dimensions(),
       QMP_get_logical_coordinates());
#else
  { int i;
    for (i = 0; i < 4; i++)
       mynode[i] = 0;
  }
#endif

  if (QOP_MDWF_init(&mdwf_state,
		    mylattice, mynetwork, mynode, primary,
		    getsub, NULL)) {
    zprint("MDWF_init() failed");
    goto end;
  }

  zprint("MDWF_init() done");

  if (QOP_MDWF_set_generic(&mdwf_params, mdwf_state, b5, c5, 0.123, 0.05)) {
    zprint("MDW_set_generic() failed");
    goto end;
  }
  zprint("MDWF_set_generic() done");

  if (do_run(mdwf_state, mdwf_params)) {
    zprint("float test failed");
    goto end;
  }

  QOP_MDWF_fini(&mdwf_state);

  zprint("Heater test finished");
  status = 0;
 end:
  QMP_finalize_msg_passing();
  return status;
}
