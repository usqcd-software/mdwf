#ifndef MARK_ADA7147B_644F_4A44_B80B_35994BC51EAE
#define MARK_ADA7147B_644F_4A44_B80B_35994BC51EAE

struct Fermion *new_fermion(int size, int ls);
struct vFermion *new_vfermion(int size, int ls, int dim);

void mk_fermion(struct Fermion *r, int size, int ls, const int *p, int off);
void show_fermion(const char *n, const struct Fermion *f, int size, int ls);

void construct_f(int esize, int ls, struct Fermion *f, double m);
void construct_vf(int esize, int ls, int width,
                  struct vFermion *dst,
                  struct Fermion *t, double m);
void construct_d(int len, double *d, double m);
#endif /* !defined(MARK_ADA7147B_644F_4A44_B80B_35994BC51EAE) */
