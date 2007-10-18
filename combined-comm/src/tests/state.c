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

void
q(set_down_pack)(struct down_pack *dp, int p, int f)
{
  xprint("%p.dp[%5d]=%5d", dp, p, f);
}

void
q(set_up_pack)(struct up_pack *dp, int p, int f, int u)
{
  xprint("%p.up[%5d]=%5d %5d", dp, p, f, u);
}

void
q(set_neighbor)(struct neighbor *np, int p, int m,
		const int f_u[Q(DIM)], int u_u,
		const int f_d[Q(DIM)], const int u_d[Q(DIM)])
{
  xprint("%p.nb[%5d]= %02x u %5d %5d %5d %5d . %5d"
	 " d %5d %5d %5d %5d . %5d %5d %5d %5d",
	 np, p, m,
	 f_u[0], f_u[1], f_u[2], f_u[3],
	 u_u,
	 f_d[0], f_d[1], f_d[2], f_d[3],
	 u_d[0], u_d[1], u_d[2], u_d[3]);
}

void
q(fix_neighbor_f_up)(struct neighbor *np, int p, int f_up, int d)
{
  xprint("%p.nb[%5d]= fix u[%d] %5d",
	 np, p, d, f_up);
}

void
q(fix_neighbor_f_down)(struct neighbor *np, int p, int f_up, int d)
{
  xprint("%p.nb[%5d]= fix d[%d] %5d",
	 np, p, d, f_up);
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
dump_eo(const char *name, struct QOP_MDWF_State *state, struct eo_lattice *eo)
{
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
}

static void
dump_state(struct QOP_MDWF_State *state)
{
  xprint("MDWF state at %p", state);
  xprint("allocate=%d", state->allocated);
  dump_eo("even", state, &state->even);
  dump_eo("odd", state, &state->odd);
}

int
main(int argc, char *argv[])
{
  struct QOP_MDWF_State *mdwf_state = NULL;
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
    zprint("arg[%d]=%s", i, argv[i]);
  if (argc != 10) {
    zprint("10 arguments expected, found %d", argc);
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
    zprint("declare_logical_top failed");
    goto end;
  }

  getv(mynode, 0, QMP_get_logical_number_of_dimensions(),
       QMP_get_logical_coordinates());

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
