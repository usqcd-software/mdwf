#include <qmp.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mdwf.h>

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
xprint(char *fmt, ...)
{
  char buffer[14096];
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
  char buffer[14096];
  va_list va;
  if (primary == 0)
    return;

  va_start(va, fmt);
  vsprintf(buffer, fmt, va);
  va_end(va);
  printf("zzz %s\n", buffer);
  fflush(stdout);
}

static void
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

static void
show_4d(const char *name, const char *part,
	struct QOP_MDWF_State *state,
	int v[])
{
  xprint("%p.%s.%s %5d %5d %5d %5d",
	 state, name, part, v[0], v[1], v[2], v[3]);
}

static void
dump_down_send(const char *name, int d,
	       struct QOP_MDWF_State *state,
	       struct eo_lattice *oe,
	       int send_size,
	       struct down_pack *send_down)
{
  int i;

  if (send_size == 0)
    return;

  for (i = 0; i < send_size; i++) {
    int f;
    q(get_down_pack)(&f, send_down, i);
    xprint("snd.%s.d[%d][%5d] = %5d", name, d, i, f);
  }
}

static void
dump_up_send(const char *name, int d,
	     struct QOP_MDWF_State *state,
	     struct eo_lattice *oe,
	     int send_size,
	     struct up_pack *send_up)
{
  int i;

  if (send_size == 0)
    return;

  for (i = 0; i < send_size; i++) {
    int f, u;
    q(get_up_pack)(&f, &u, send_up, i);
    xprint("snd.%s.u[%d][%5d] = %5d u %5d", name, d, i, f, u);
  }
}

static void
dump_eo(const char *name, struct QOP_MDWF_State *state,
	struct eo_lattice *eo,
	struct eo_lattice *oe)
{
  int i, d;

  xprint("%p.%s face %6d", state, name, eo->face_size);
  xprint("%p.%s body %6d", state, name, eo->body_size);
  xprint("%p.%s full %6d", state, name, eo->full_size);
  xprint("%p.%s neighbor %p", state, name, eo->body_neighbor);
  show_4d(name, "local.lo", state, eo->local->lo);
  show_4d(name, "local.hi", state, eo->local->hi);
  show_4d(name, "local.dx", state, eo->local->dx);
  show_4d(name, "down_snd", state, eo->send_down_size);
  show_4d(name, "down_rcv", state, eo->receive_down_size);
  show_4d(name, "up_snd", state, eo->send_up_size);
  show_4d(name, "up_rcv", state, eo->receive_up_size);
  for (i = 0; i < eo->full_size; i++) {
    int m, u_u;
    int x[4];
    int f_u[4], f_d[4], u_d[4];
    q(get_neighbor)(&m, f_u, &u_u, f_d, u_d, eo->body_neighbor, i);
    q(l2v)(x, eo->local, eo->lx2v[i]);
    xprint("%s.nb %5d {%5d, %5d, %5d, %5d}: %02x"
	   " f. %5d %5d %5d %5d  :: %5d,"
	   " b. %5d %5d %5d %5d  :: %5d %5d %5d %5d",
	   name, i,
	   x[0], x[1], x[2], x[3],
	   m,
	   f_u[0], f_u[1], f_u[2], f_u[3], u_u,
	   f_d[0], f_d[1], f_d[2], f_d[3],
	   u_d[0], u_d[1], u_d[2], u_d[3]);
  }
  for (d = 0; d < 4; d++) {
    dump_down_send(name, d, state, oe,
		   eo->send_down_size[d],
		   eo->down_pack[d]);
    dump_up_send(name, d, state, oe,
		 eo->send_up_size[d],
		 eo->up_pack[d]);
  }
}

static void
dump_state(struct QOP_MDWF_State *state)
{
  int i;

  xprint("MDWF state at %p", state);
  xprint("allocate = %d", state->allocated);
  xprint("v2lx     = %p", state->v2lx);
  xprint("volume   = %d", state->volume);
  for (i = 0; i < state->volume; i++) {
    int x[4];
    q(l2v)(x, &state->local, state->lx2v[i]);
    xprint(" u[%5d]: { %5d, %5d, %5d, %5d}",
	   i, x[0], x[1], x[2], x[3]);
  }
  dump_eo("even", state, &state->even, &state->odd);
  dump_eo("odd", state, &state->odd, &state->even);
}

int
main(int argc, char *argv[])
{
  struct QOP_MDWF_State *mdwf_state = NULL;
  QMP_thread_level_t qt = QMP_THREAD_SINGLE;
  int status = 1;
  int vx[4];
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
    zprint("arg[%d]=%s", i, argv[i]);
  if (argc != 14) {
    zprint("14 arguments expected, found %d", argc);
    QMP_finalize_msg_passing();
    return 1;
  }

  for (i = 0; i < 4; i++) {
    mynetwork[i] = atoi(argv[i+1]);
    mylocal[i] = atoi(argv[i+5]);
    vx[i] = atoi(argv[i+10]);
    mylattice[i] = mylocal[i] * mynetwork[i];
  }
  mylocal[4] = mylattice[4] = atoi(argv[9]);

  zshowv4("network", mynetwork);
  zshowv5("local lattice", mylocal);
  zshowv5("lattice", mylattice);

  getv(mynode, 0, 4, vx);

  xshowv("node", mynode);

  if (QOP_MDWF_init(&mdwf_state,
		    mylattice, mynetwork, mynode, primary,
		    getsub, NULL)) {
    zprint("MDWF_init() failed");
    goto end;
  }

  zprint("MDWF_init() done");

  dump_state(mdwf_state);

  QOP_MDWF_fini(&mdwf_state);

  status = 0;
 end:
  QMP_finalize_msg_passing();
  return status;
}
