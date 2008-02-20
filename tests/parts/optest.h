#ifndef _optest_h
#define _optest_h
#define NELEM(x) (sizeof (x) / sizeof ((x)[0]))

extern int fermion_pos[5];
extern int fermion_color;
extern int fermion_dirac;
extern int fermion_reim;

extern char *op_name;
extern int operator(void);

extern double read_gauge(int dir,
			 const int pos[4],
			 int a, int b,
			 int re_im,
			 void *env);
extern double read_fermion(const int pos[5],
			   int c, int d,
			   int re_im,
			   void *env);
extern void write_fermion(const int pos[5],
			  int c, int d,
			  int re_im,
			  double value,
			  void *env);
extern void setup_bc(void);

extern struct QOP_MDWF_State *state;
extern struct QOP_MDWF_Parameters *params;
extern struct QOP_MDWF_Fermion *result;
extern struct QOP_MDWF_Fermion *fermion;
extern struct QOP_MDWF_Gauge *gauge;

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
