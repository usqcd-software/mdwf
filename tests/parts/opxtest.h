#ifndef _opxtest_h
#define _opxtest_h
#define NELEM(x) (sizeof (x) / sizeof ((x)[0]))

unsigned int sum_init(int v);
unsigned int sum_add(unsigned int state, int v);
double sum_fini(unsigned int state);

extern char *op_a_name;
extern char *op_b_name;
extern int operator_a(void);
extern int operator_b(void);

extern double read_gauge(int dir,
			 const int pos[4],
			 int a, int b,
			 int re_im,
			 void *env);

extern double read_fermion_a(const int pos[5],
			     int c, int d,
			     int re_im,
			     void *env);
extern double read_fermion_b(const int pos[5],
			     int c, int d,
			     int re_im,
			     void *env);

extern struct QOP_MDWF_State *state;
extern struct QOP_MDWF_Parameters *params;
extern struct QOP_MDWF_Fermion *fermion_a;
extern struct QOP_MDWF_Fermion *fermion_b;
extern struct QOP_MDWF_Gauge *gauge;
extern unsigned int seed_u;
extern unsigned int seed_a;
extern unsigned int seed_b;

extern int self;
extern int primary;
extern int lattice[5];
extern int network[4];
extern int node[4];
extern double b5[128];
extern double c5[128];

void zflush(void);
void xprint(char *fmt, ...);
void zprint(char *fmt, ...);

#endif
