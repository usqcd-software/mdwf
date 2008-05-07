#ifndef SAMPLE_COMMON_H
#define SAMPLE_COMMON_H
extern double U_scale;
extern int lattice[5];
extern int network[4];
extern int this_node[4];
extern int primary_p;

extern double M;
extern double m_5;
extern double kappa;

extern int max_iterations;
extern double min_epsilon;
extern unsigned U_seed;
extern unsigned rhs_seed;
extern unsigned sol_seed;
extern unsigned options;

double read_gauge(int d, const int p[4], int a, int b, int re_im, void *env);
double read_fermion(const int p[5], int c, int d, int re_im, void *env);
void write_fermion(const int p[5], int c, int d, int re_im, double v, void *e);
int init_qmp(int argc, char *argv[], const char *name, char prec,
	     int *count, double **shift);
void fini_qmp(void);
void zprint(const char *who, const char *fmt, ...);
void get_sublattice(int lo[4], int hi[4], const int node[4], void *env);
void report_performance(struct QOP_MDWF_State *state, char *name);
void report_time(struct QOP_MDWF_State *state, char *name);

#endif
