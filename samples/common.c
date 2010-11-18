/* This is a file of helper functions and data for sample solvers.
 */
#include <qmp.h>
#include <qop-mdwf3.h>
#include <stdarg.h>
#include <stdlib.h>
#include "common.h"

static int qmp_inited = 0;

double U_scale;
int lattice[5] = { -1, -2, -3, -4, -5}; /* fly tape begins */
int network[4] = { -6, -7, -8, -9};
int this_node[4] = { -10, -11, -12, -13};
int neighbor_up[4];
int neighbor_down[4];
/*
int neighbor_up[4] = { -1, -2, -3, -4 };
int neighbor_down[4] = { -1, -2, -3, -4 };
*/
 /* fly tape ends */
int primary_p;

double M;
double m_5;
double kappa;

int max_iterations;
double min_epsilon;
unsigned U_seed;
unsigned rhs_seed;
unsigned sol_seed;
unsigned options;

static int
node_4to1(const int n4[])
{
	int i, n;
	for (n = 0, i = 4; i--;) {
		n = n * network[i] + n4[i];
	}
	return n;
}

static void
node_1to4(int n4[], int node)
{
	int i;
	for (i = 0; i < 4; i++) {
		n4[i] = node % network[i];
		node = node / network[i];
	}
}

/* Not very good random numbers */
#define NELEM(x) (sizeof (x) / sizeof (x[0]))

static unsigned int const table[] = {
    0x575c3abc, 0x66a6de09, 0x86575ad4, 0x34bcd1e5,
    0x50f3387f, 0x6352990d, 0x93532d54, 0x3848eb41,
    0x93372ae2, 0xc7ee6706, 0xde90a3e8, 0x8637a642,
    0x4faa1a42, 0x531c15ed, 0xb016512b, 0xdfb23921,
    0x37b2bf1a, 0xee91941f, 0x1395ef23, 0x302481a4,
    0xd5e714de, 0x177cf612, 0xa75b371d, 0x9b87e32c,
    0xc6260af1, 0x88f561f1, 0x9a943091, 0x8ca94c6c,
    0x9cb5ecb3, 0x3b99aac2, 0x1cda8a11, 0xcbdf2be9,
    0x378e6a21, 0xbf9c60cd, 0x78fcf142, 0x0ffdac53,
    0x98cd5852, 0xe4a55aec, 0x4bc6be5a, 0x9f735c12,
    0x3a33ad62, 0xd03d0eb2, 0xcff3aeab, 0x30edf951,
    0xff4d23f1, 0xbdd81b0a, 0x785b5e2e, 0x592de686,
    0x9a8b2865, 0x18a9ca91, 0x526315d9, 0x6137d1e8,
    0x990134cd, 0x29d36be8, 0x62a49721, 0x53b6af9b,
    0x867a9ea1, 0x3dd656fb, 0x979ff13b, 0xbb2d9d1f,
    0x14fc90da, 0x590acdaa, 0xd785131a, 0xbb6900d4,
    0xd4b2f271, 0xe240168b, 0x1dd30da7, 0x40a98d33,
    0x9aa26a00, 0x532c6f13, 0x22ffe346, 0xcd2e49ca,
    0x39d1ecb4, 0x33572670, 0xdbbfec86, 0x729e2f73,
    0x5e35cb1e, 0x68a70192, 0xc8d6f576, 0x30e965ee,
    0x930679b7, 0xd7e6c2dd, 0xc67477b5, 0x7a6b17e0,
    0x7e0acfaa, 0x532d7403, 0x338dc81a, 0x14be5463,
    0x3ba9d6da, 0x4a5962a7, 0xc8778225, 0x3df3dd6c,
    0x5adbd365, 0x04ff3fcd, 0x74094302, 0xef8f7724,
    0xb57107c1, 0xda9a5e10, 0x686e496c, 0xa39eb159,
    0xd42bf574, 0x47fd5be1, 0xaa102ed7, 0xae6be4ce,
    0x80bb8eba, 0x7940b332, 0x77a9f951, 0x66a2e3ee,
    0x6e8574af, 0x0db6f940, 0xea9556d5, 0x941667c1,
    0xc329548f, 0xd5cfa91a, 0x0b74735c, 0x6aadc0c0,
    0xbb0c11ea, 0x65feb3eb, 0xcc8c19a3, 0xb1f6d56a,
    0x4307ef7b, 0x0a9c92ad, 0xcf2d81e3, 0xccd24a49,
    0xb4db001b, 0x87c2219b, 0x91d0d127, 0x938cb6d6,
    0xcca896d4, 0x47b6e4f8, 0x76408a35, 0xe27b2efe,
    0x1e4a72e3, 0x7f458dd0, 0x57076fa2, 0x018ee87b,
    0x1664f1fd, 0x793343b4, 0xb283a17b, 0x2e369ce8,
    0xfe7cacf2, 0x8af45cc1, 0xfa071d01, 0x21c3511f,
    0xdaa5702e, 0x641ef2cc, 0x198e97f7, 0x5ffbd66b,
    0xa0a65d77, 0x951f0b53, 0x6806fcf5, 0xf6cc841e,
    0x5ca947e3, 0x23aa2651, 0x28c7394f, 0xd228d400,
    0xa277be27, 0x5c355284, 0x467fd143, 0xa5f0201e,
    0x1e050f5f, 0xc0d9f3c2, 0x5504f293, 0x1991f2d6,
    0x895fc215, 0xd255d0f9, 0xcc5acd5e, 0x01ff4e98,
    0xfc30bec8, 0x1080c9a2, 0x8854cc32, 0xe8adfa26,
    0x1bbdc1c5, 0x7a7020f4, 0x393c6848, 0xe3132e2a,
    0xd7b7853d, 0x0b0bcb24, 0x95ceaffb, 0x60341a5a,
    0xe036c7bd, 0x928e684a, 0x9b6a4929, 0xe12ac3c2,
    0x16a7f7a7, 0x0891e1b2, 0xc93f64e6, 0x0c286347,
    0xe860182f, 0xf5946550, 0x3ba46dda, 0x54ae7a82,
    0xc6b4f3ee, 0xee7dc63b, 0xdcc7f34e, 0x805b20d2,
    0xb13dcd86, 0xbd1799ae, 0x7b472dc1, 0x11c148f2,
    0x2a7ef148, 0x17fa1c47, 0x1dacdf5a, 0x78deef13,
    0x516aa8ef, 0xfa43f1ae, 0xc7f1bd2f, 0xb7b08044,
    0xdfae280e, 0x4bf43476, 0x99a5875e, 0x43737c43,
    0x6f408478, 0x62cd4cfc, 0x860f8c29, 0x70e766bf,
    0x188e67a9, 0xe02905a4, 0x432c9a27, 0xca7d64ce,
    0x2e25b4ad, 0x7e3a1908, 0x9cd8549d, 0x002266e3,
    0x98ba14c2, 0x3f18dbc6, 0x6b9eb0c8, 0xd46f97ad,
    0x52ad030a, 0xcaf5b07f, 0x471bdf3b, 0xe2bc4855,
    0x57464d8b, 0x7ba1e756, 0xff9ce231, 0xca850d81,
    0xe6ed2cd5, 0x46efaaad, 0x643e3dd6, 0x7ccc48a3,
    0x4991072c, 0x2c602e37, 0xb84f6af4, 0x1ece9dae,
    0xca89f5a1, 0x9b9e72bb, 0x50b4742b, 0x000ec926,
    0x57202775, 0xf4f9a05e, 0x4b0f578b, 0xd5fda952,
    0x284b7d80, 0x2872b3fa, 0xbeb9553b, 0x8cd00b67,
    0xb7007421, 0xff532635, 0x76eb8d5e, 0x0cef9e8c,
    0x4665bb19, 0xf565675c, 0xbfd4f905, 0x9629ca3c,
    0x53d86547, 0x837869d8, 0x967059d1, 0x901d5ed7,
    0x1f9b40ce, 0x800c28a7, 0xb38ee9b9, 0x9d2a169e,
    0x353bf221, 0xe8dc1c44, 0x16894de5, 0xaaa35e87,
    0x81e4cc91, 0xf7621e3b, 0x70c4e04b, 0x55dc261d,
    0x70f86b7b, 0x191622f5, 0xc2d85901, 0xceecefb4,
    0x51aa671c, 0x3ff5f991, 0x53f877b0, 0x76c862f5,
    0xa053e4c5, 0x61c4e4a6, 0x451389ba, 0x52fa9216,
    0x21179e87, 0xc3c3e17b, 0x2e26d1dc, 0x6e0a188a,
    0xe2270065, 0x04bce1c4, 0xcf337da7, 0x1ae7b67c,
    0x27f2383e, 0xbc2ad6cd, 0x5d319aa2, 0xd48df243,
    0x968cedbd, 0x562bfb88, 0x6852f2a5, 0x85c101ec,
    0xd8b4f830, 0x88715750, 0x46f094b8, 0x49e7329b,
    0x6a4a5a0c, 0x46790032, 0x5c9e8162, 0xca1792a4,
    0xb72739f0, 0xfcdf1da4, 0x92300afc, 0xb0b948ae,
    0xa91dc6da, 0xaf5722f1, 0xe5a3625a, 0x0c0305a9,
    0xf6bf6316, 0xe2bdf118, 0xc552194b, 0x44b6c120,
    0x87f56de0, 0x66719b06, 0xd939d161, 0x9c16383b,
    0x7195fd40, 0x5cfc6c55, 0xe9c0afe4, 0xb4d70e79,
    0x0ebe311b, 0x744fd740, 0x8da9bad1, 0xb81f5f47,
    0x1e288180, 0xe9af2f19, 0x4dee6dd5, 0x40c1a722,
    0x7b579e93, 0x5a4e3a1d, 0x372eef27, 0x6aaddb1d,
    0x0ae8386e, 0x77c89b43, 0x363b183c, 0x733c560a,
    0xf2eb0d71, 0x283c895e, 0xda3698bc, 0x86565107,
    0xf9c51c8b, 0x12d35110, 0xb93f3cce, 0x38501055,
    0x1842849a, 0x7f68c937, 0x7c21df77, 0x40ed1a3b,
    0x2ae5a413, 0xa04385f7, 0x4ee9c482, 0x58b34205,
    0x806eb0bf, 0x057f828d, 0xb9d44d37, 0xfe748fa9,
    0xb38dc19e, 0x8119d698, 0x50c57afa, 0x24104f2b,
    0x600dd625, 0xef20409e, 0x8c0daea8, 0xfa4cf5cc,
    0xc9310ec7, 0xba16dcf2, 0x412d6c8e, 0x55914283,
    0xee74714e, 0x3c95e47c, 0x9df2c23a, 0x738504f1,
    0xe5fa4256, 0xee70f5d8, 0x1729f855, 0x91fe7113,
    0x04e0b7f1, 0x263b11a8, 0x48e705d9, 0x94647bdf,
    0x463e0c58, 0x841384e6, 0x9b01a19d, 0xa370260a,
    0x88b6a407, 0x2bbd0c49, 0x92239b04, 0x70ae3b56,
    0x5d3ccf56, 0x8147e442, 0xcbccaf1e, 0x7274d6fc};

unsigned int
sum_init(int v)
{
    return v;
}

unsigned int
sum_add(unsigned int state, int v)
{
    return state ^ table[((state >> (v & 7)) + v + 1) % NELEM(table)];
}

double
sum_fini(unsigned int state)
{
    return state * 4.6566128730773925781e-10 - 1;
}

/* get local sublattice low and high bounds */
void
get_sublattice(int lo[4], int hi[4], int node, void *env)
{
    int i;
	int n4[4];

	node_1to4(n4, node);
    for (i = 0; i < 4; i++) {
	lo[i] = (lattice[i] * n4[i]) / network[i];
	hi[i] = (lattice[i] * (n4[i] + 1)) / network[i];
    }
}

/* empty fermion writer */
void
write_fermion(const int pos[5],
	      int c, int d,
	      double v_re,
	      double v_im,
	      void *env)
{
}

/* read fermion. Produce a random number from -1 to +1 for each component */
void
read_fermion(double *v_re,
			 double *v_im,
			 const int pos[5],
			 int c, int d,
			 void *env)
{
    int i;
    unsigned int v = sum_init(*(int *)env);

    v = sum_add(v, c);
    v = sum_add(v, d);
    for (i = 0; i < 5; i++) {
		v = sum_add(v, pos[i]);
		v = sum_add(v, *(int *)env);
    }
	*v_re = sum_fini(v);
    v = sum_add(v, 1);
	*v_im = sum_fini(v);
}

/* read gauge. Produce a random matrix plus a unit */
void
read_gauge(double *v_re,
		   double *v_im,
		   int dir,
		   const int pos[4],
		   int a, int b,
		   void *env)
{
    int i;
    unsigned int v = sum_init(*(int *)env);
    double x;

    v = sum_add(v, a);
    v = sum_add(v, b);
    for (i = 0; i < 4; i++) {
	v = sum_add(v, pos[i]);
	v = sum_add(v, dir);
    }
    v = sum_add(v, *(int *)env);
    x = U_scale * sum_fini(v);
    if (a == b)
		x = x + 1;
	*v_re = x;
    v = sum_add(v, 1);
	*v_im = sum_fini(v);
}

/* print stuff on primary node only */
void
zprint(const char *who, const char *fmt, ...)
{
    char buffer[4096];
    va_list va;

    if (primary_p == 0)
	return;
    
    va_start(va, fmt);
    vsnprintf(buffer, sizeof (buffer) - 1, fmt, va);
    va_end(va);

    printf("%s: %s\n", who, buffer);
}

int
init_qmp(int argc, char *argv[], const char *who, char prec,
	 int *count, double **shift)
{
    QMP_thread_level_t qt = QMP_THREAD_SINGLE;
    const int *qmp_addr;
    const char *marker;
    int qmp_dim;
    int i;

    if (QMP_init_msg_passing(&argc, &argv, qt, &qt) != QMP_SUCCESS) {
	fprintf(stderr, "QMP_init() failed\n");
	return 1;
    }
    primary_p = QMP_is_primary_node();
    
    zprint(who, "Starting precision (%c) example", prec);
    zprint(who, "MDWF version: %s", QOP_MDWF_version());

    /* get parameters */
    if (((count == NULL) && (argc != 21)) ||
	((count != NULL) && (argc <= 22))) {
	zprint(who,
	       "usage: solver-mxm Nx Ny Nz Nt Lx Ly Lz Lt Ls"
	       " M m_5 kappa"
	       " guess-seed source-seed gauge-seed gauge-scale"
	       " max-iterations min-epsilon verbose-p marker%s",
	       count? " shift ...": "");
	goto err_0;
    }
    for (i = 0; i < 4; i++)
	network[i] = atoi(argv[i+1]);
    for (i = 0; i < 5; i++)
	lattice[i] = atoi(argv[i+5]);
    M = atof(argv[10]);
    m_5 = atof(argv[11]);
    kappa = atof(argv[12]);
    sol_seed = atoi(argv[13]);
    rhs_seed = atoi(argv[14]);
    U_seed = atoi(argv[15]);
    U_scale = atof(argv[16]);
    max_iterations = atoi(argv[17]);
    min_epsilon = atof(argv[18]);
    options = atoi(argv[19]);
    marker = argv[20];
    if (options)
	options = QOP_MDWF_LOG_EVERYTHING;
    if (count) {
	int i;

	*count = argc - 21;
	*shift = malloc(*count * sizeof (double));
	if (*shift == 0) {
	    zprint(who, "not enough memory");
	    return 1;
	}
	for (i = 0; i < *count; i++)
	    (*shift)[i] = atof(argv[21 + i]);
    }


    /* print what we've got */
    zprint(who, "QMP nodes: %d", QMP_get_number_of_nodes());
    zprint(who, "marker: %s", marker);
    zprint(who, "requesting %d nodes",
	   network[0] * network[1] * network[2] * network[3]);
    zprint(who, "network: %d %d %d %d",
	   network[0], network[1], network[2], network[3]);
    zprint(who, "lattice: %d %d %d %d %d",
	   lattice[0], lattice[1], lattice[2], lattice[3], lattice[4]);
    zprint(who, "M   = %8.6f", M);
    zprint(who, "m_5 = %8.6f", m_5);
    zprint(who, "kappa = %8.6f", kappa);
    zprint(who, "seed A = %d", sol_seed);
    zprint(who, "seed B = %d", rhs_seed);
    zprint(who, "seed C = %d", U_seed);
    zprint(who, "scale U = %g", U_scale);
    zprint(who, "max iterations = %d", max_iterations);
    zprint(who, "min epsilon = %e", min_epsilon);
    zprint(who, "options = %d", options);
    zprint(who, "count=%p, shift=%p", count, shift);
    if (count) {
	int i;

	zprint(who, "shift count = %d", *count);
	for (i = 0; i < *count; i++) {
	    zprint(who, "shift[%d] = %8.6f", i, (*shift)[i]);
	}
    }

    if (QMP_declare_logical_topology(network, 4) != QMP_SUCCESS) {
	zprint("ERROR", "declare_logical_topology failed");
	goto err_0;
    }

    qmp_dim = QMP_get_logical_number_of_dimensions();
    qmp_addr = QMP_get_logical_coordinates();
	{
		/* setup neighbors */
		int d, i, n4[4], v[4];
		node_1to4(n4, QMP_get_node_number());
		for (d = 0; d < 4; d++) {
			for (i = 0; i < 4; i++)
				v[i] = n4[i];
			if (v[d] + 1 == network[d])
				v[d] = 0;
			else
				v[d] = v[d] + 1;
			neighbor_up[d] = node_4to1(v);
			for (i = 0; i < 4; i++)
				v[i] = n4[i];
			if (v[d] == 0)
				v[d] = network[d] - 1;
			else
				v[d] = v[d] - 1;
			neighbor_down[d] = node_4to1(v);
		}
	}

    qmp_inited = 1;
    return 0;
err_0:
    QMP_finalize_msg_passing();
    return 1;
}


void
fini_qmp(void)
{
    if (qmp_inited == 0)
	return;

    QMP_finalize_msg_passing();
}


void
report_performance(struct QOP_MDWF_State *state, char *name)
{
    double total_time;
    long long total_flops;

    /* Get statistics from the MDWF layer */
    if (QOP_MDWF_performance(&total_time,
			     &total_flops,
			     NULL, NULL,
			     state)) {
	zprint(name, "perf() returned error? This is strange but possible");
    } else if (total_time == 0.0) {
	zprint(name, "too short time interval");
    } else {
	zprint(name, "total time %.3f sec, performance %.3f MFlop/s",
	       total_time, total_flops * 1e-6 / total_time);
    }
}

void
report_time(struct QOP_MDWF_State *state, char *name)
{
    double total_time;

    /* Get statistics from the MDWF layer */
    if (QOP_MDWF_performance(&total_time,
			     NULL, NULL, NULL,
			     state)) {
	zprint(name, "perf() returned error? This is strange but possible");
	return;
    }
    zprint(name, "total time %.3f sec", total_time);
}
