#include <qmp.h>
#include <stdlib.h>
#include <stdarg.h>

#define NELEM(x) (sizeof (x) / sizeof (x[0]))

static int self;
static int primary;
static int mylattice[5];
static int mylocal[5];
static int mynetwork[4];
static int mynode[4];
static double b5[128];
static double c5[128];

void
zflush(void)
{
    if (primary == 0)
	return;
    fflush(stdout);
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
getsub(int lo[4], int hi[4], const int node[4], void *env)
{
    int i;

    for (i = 0; i < 4; i++) {
	lo[i] = (mylattice[i] * node[i]) / mynetwork[i];
	hi[i] = (mylattice[i] * (node[i] + 1)) / mynetwork[i];
    }
}

static double
read_gauge(int dir,
	   const int pos[4],
	   int a,
	   int b,
	   int re_im,
	   void *env)
{
    if (a == b && re_im == 0)
	return 1.0;
    else
	return 0.0;
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

struct QOP_MDWF_State;
static int report_perf(const char *name, struct QOP_MDWF_State *state,
		       int iter, double epsilon);

#define QOP_MDWF_DEFAULT_PRECISION 'F'
#define do_xxx do_f3
#include "do-cg.c"

#define QOP_MDWF_DEFAULT_PRECISION 'D'
#define do_xxx do_d3
#include "do-cg.c"

static int
report_perf(const char *name, struct QOP_MDWF_State *state,
	    int iter, double epsilon)
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
    rounds ++;
    
    if (t_sec < 10)
	return 1;
    
    zprint("CG perf(%s @ %d): %g MF/sec, snd %g MB/sec, rcv %g MB/sec "
	   "[%g %lld %lld %lld] %d iterations, espilon %g",
	   name, rounds,
	   1e-6 * t_flops/t_sec,
	   1e-6 * t_send/t_sec,
	   1e-6 * t_receive/t_sec,
	   t_sec, t_flops, t_send, t_receive,
	   iter, epsilon);
    zflush();
    return 0;
}

static void
get_vector(int v[4], int def, int dim, const int d[])
{
    int i;
    
    for (i = 0; i < 4; i++)
	v[i] = def;
    for (i = 0; i < dim && i < 4; i++)
	v[i] = d[i];
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
    if (argc != 10) {
	zprint("10 arguments expected, found %d", argc);
	zprint("Usage: conjgrad Nx Ny Nz Nt lx ly lz lt Ls");
	QMP_finalize_msg_passing();
	return 1;
    }

    for (i = 0; i < 4; i++) {
	mynetwork[i] = atoi(argv[i+1]);
	mylocal[i] = atoi(argv[i+5]);
	mylattice[i] = mylocal[i] * mynetwork[i];
    }
    mylocal[4] = mylattice[4] = atoi(argv[9]);
    
    zprint("lattice %d %d %d %d %d",
	   mylattice[0],
	   mylattice[1],
	   mylattice[2],
	   mylattice[3],
	   mylattice[4]);
    zprint("network %d %d %d %d",
	   mynetwork[0],
	   mynetwork[1],
	   mynetwork[2],
	   mynetwork[3]);
    
    if (QMP_declare_logical_topology(mynetwork, 4) != QMP_SUCCESS) {
	zprint("declare_logical_top failed");
	goto end;
    }

    get_vector(mynode, 0,
	       QMP_get_logical_number_of_dimensions(),
	       QMP_get_logical_coordinates());

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

    if (do_f3(mdwf_state, mdwf_params, "cg_f3()")) {
	zprint("float test failed");
	goto end;
    }
    if (do_d3(mdwf_state, mdwf_params, "cg_d3()")) {
	zprint("double test failed");
	goto end;
    }

    QOP_MDWF_fini(&mdwf_state);
    zprint("Conjugate gradient test finished");
    status = 0;
 end:
    QMP_finalize_msg_passing();
    return status;
}
