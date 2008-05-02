#define QOP_MDWF_DEFAULT_PRECISION 'F'
#include <qop-mdwf3.h>
#include <mdwf.h>
#include <qmp.h>
#include <stdlib.h>
#include <stdarg.h>

#define NELEM(x) (sizeof (x) / sizeof (x[0]))

static int self;
static int primary;
static unsigned int sU;
static unsigned int sX;
static unsigned int sY;
static double M;
static double m_5;
static int mylattice[5];
static int mynetwork[4];
static double b5[128];
static double c5[128];
static int mynode[4];

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
setup_bc(unsigned int sB, unsigned int sC)
{
    int i;
    unsigned int v;
    
    for (i = 0; i < mylattice[4]; i++) {
	v = sum_init(sB);
	v = sum_add(v, i);
	v = sum_add(v, sC);
	b5[i] = sum_fini(v);
	v = sum_init(sC);
	v = sum_add(v, i);
	v = sum_add(v, sB);
	v = sum_add(v, sC);
	c5[i] = sum_fini(v);
    }
}

void
zflush(void)
{
    if (primary == 0)
	return;
    fflush(stdout);
}

void
zzprint(char *fmt, ...)
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

#include "do-cg.c"

static double
read_gauge(int dir,
	   const int pos[4],
	   int a, int b,
	   int re_im,
	   void *env)
{
    int i;
    unsigned int v = sum_init(*(int *)env);
    v = sum_add(v, re_im);
    v = sum_add(v, a);
    v = sum_add(v, b);
    for (i = 0; i < 4; i++) {
	v = sum_add(v, pos[i]);
	v = sum_add(v, dir);
    }
    v = sum_add(v, *(int *)env);
    return sum_fini(v);
}

static double
read_fermion(const int pos[5],
	     int c, int d,
	     int re_im,
	     void *env)
{
    int i;
    unsigned int v = sum_init(*(int *)env);

    v = sum_add(v, c);
    v = sum_add(v, d);
    v = sum_add(v, re_im);
    for (i = 0; i < 5; i++) {
	v = sum_add(v, pos[i]);
	v = sum_add(v, *(int *)env);
    }
    return sum_fini(v);
}


static int
run_it(struct QOP_MDWF_State *state, struct QOP_MDWF_Parameters *params,
       int max_iter, double eps)
{
    char *name = "run_it";
    struct QOP_MDWF_Gauge *U;
    struct QOP_MDWF_Fermion *rhs;
    struct QOP_MDWF_Fermion *guess;
    struct QOP_MDWF_Fermion *R;
    int out_iter;
    double out_eps;

    zzprint("%s: starting conjugate gradient test in precision %c",
	   name, QOP_MDWF_DEFAULT_PRECISION);

    if (QOP_MDWF_import_gauge(&U, state, read_gauge, &sU)) {
	zzprint("%s: import gauge failed", name);
	goto no_U;
    }
    if (QOP_MDWF_import_fermion(&rhs, state, read_fermion, &sY)) {
	zzprint("%s: import rhs failed", name);
	goto no_rhs;
    }
    if (QOP_MDWF_import_fermion(&guess, state, read_fermion, &sX)) {
	zzprint("%s: import initial guess failed", name);
	goto no_guess;
    }
    if (QOP_MDWF_allocate_fermion(&R, state)) {
	zzprint("%s: allocating solution failed", name);
	goto no_R;
    }

    cg(R, &out_iter, &out_eps,
       params, guess, U, rhs, max_iter, eps,
       QOP_MDWF_LOG_EVERYTHING);

    QOP_MDWF_free_fermion(&R);
    QOP_MDWF_free_fermion(&guess);
    QOP_MDWF_free_fermion(&rhs);
    QOP_MDWF_free_gauge(&U);
    zzprint("%s: end max_iter=%d, esp=%g", name, max_iter, eps);
    return 0;

no_R:
    QOP_MDWF_free_fermion(&guess);
no_guess:
    QOP_MDWF_free_fermion(&rhs);
no_rhs:
    QOP_MDWF_free_gauge(&U);
no_U:
    return 1;
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
    unsigned int sB;
    unsigned int sC;
    int max_iter;
    double min_prec;
    int status = 1;
    int i;
    
    if (QMP_init_msg_passing(&argc, &argv, qt, &qt) != QMP_SUCCESS) {
	fprintf(stderr, "QMP_init() failed\n");
	return 1;
    }
    
    self = QMP_get_node_number();
    primary = QMP_is_primary_node();
    if (argc != 19) {
	zzprint("17 arguments expected, found %d", argc);
	zzprint("Usage: conjgrad M m sB sC "
	       "Nx Ny Nz Nt Lx Ly Lz Lt Ls sU sX sY iter eps");
	QMP_finalize_msg_passing();
	return 1;
    }

    M = atof(argv[1]);
    m_5 = atof(argv[2]);
    sB = atoi(argv[3]);
    sC = atoi(argv[4]);
    for (i = 0; i < 4; i++) {
	mynetwork[i] = atoi(argv[i+5]);
	mylattice[i] = atoi(argv[i+9]);
    }
    mylattice[4] = atoi(argv[13]);
    sU = atoi(argv[14]);
    sX = atoi(argv[15]);
    sY = atoi(argv[16]);
    max_iter = atoi(argv[17]);
    min_prec = atof(argv[18]);
    
    zzprint("lattice %d %d %d %d %d",
	   mylattice[0],
	   mylattice[1],
	   mylattice[2],
	   mylattice[3],
	   mylattice[4]);
    zzprint("network %d %d %d %d",
	   mynetwork[0],
	   mynetwork[1],
	   mynetwork[2],
	   mynetwork[3]);
    zzprint("M = %8.6f", M);
    zzprint("m = %8.6f", m_5);
    zzprint("iter = %d", max_iter);
    zzprint("epsilon = %e", min_prec);
    zzprint("sB = %d", sB);
    zzprint("sC = %d", sC);
    zzprint("sU = %d", sU);
    zzprint("sX = %d", sX);
    zzprint("sY = %d", sY);
    
    setup_bc(sB, sC);
    for (i = 0; i < mylattice[4]; i++)
	zzprint("  [%2d] b = %16.8f   c = %16.8f", i, b5[i], c5[i]);

    if (QMP_declare_logical_topology(mynetwork, 4) != QMP_SUCCESS) {
	zzprint("declare_logical_top failed");
	goto end;
    }

    get_vector(mynode, 0,
	       QMP_get_logical_number_of_dimensions(),
	       QMP_get_logical_coordinates());

    if (QOP_MDWF_init(&mdwf_state,
		      mylattice, mynetwork, mynode, primary,
		      getsub, NULL)) {
	zzprint("MDWF_init() failed");
	goto end;
    }
    zzprint("MDWF_init() done");

    if (QOP_MDWF_set_generic(&mdwf_params, mdwf_state, b5, c5, 0.123, 0.05)) {
	zzprint("MDW_set_generic() failed");
	goto end;
    }
    zzprint("MDWF_set_generic() done");

    if (run_it(mdwf_state, mdwf_params, max_iter, min_prec)) {
	zzprint("test failed");
	goto end;
    }

    QOP_MDWF_free_parameters(&mdwf_params);
    QOP_MDWF_fini(&mdwf_state);
    zzprint("Conjugate gradient test finished");
    status = 0;
 end:
    QMP_finalize_msg_passing();
    return status;
}
