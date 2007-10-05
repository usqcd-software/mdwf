#include <qmp.h>
#include <stdarg.h>
#include <qop-mdwf3.h>

static int self;
static int primary;
static int mynetwork[4];
static int mynode[4];

static
xprint(char *fmt, ...)
{
  char buffer[4096];
  va_list va;

  va_start(va, fmt);
  vsprintf(buffer, fmt, va);
  va_end(va);
  printf("[%04d] %s\n", self, buffer);
  
}

static void
zprint(char *fmt, ...)
{
  va_list va;
  if (primary == 0)
    return;

  va_start(va, fmt);
  vprintf(fmt, va);
  va_end(va);
}

static void
show_coor(const char *name, int d, const int *dim, const int *coord)
{
  int i;
  int vd[4];
  int vc[4];

  for (i = 0; i < 4; i++)
    vd[i] = 1, vc[i] = 0;
  for (i = 0; i < d; i++)
    vd[i] = dim[i], vc[i] = coord[i];
      
  xprint("%s network: %d %d %d %d", name, vd[0], vd[1], vd[2], vd[3]);
  xprint("%s node: %d %d %d %d", name, vc[0], vc[1], vc[2], vc[3]);
}

int
main(int argc, char *argv[])
{
  QMP_thread_level_t qt = QMP_THREAD_SINGLE;
  int net_dim;
  const int *network;
  const int *node;
  int i;

  if (QMP_init_msg_passing(&argc, &argv, qt, &qt) != QMP_SUCCESS) {
    fprintf(stderr, "QMP_init() failed\n");
    return 1;
  }

  self = QMP_get_node_number();
  primary = QMP_is_primary_node();
  for (i = 0; i < argc; i++)
    zprint("arg[%d]=%s\n", i, argv[i]);

  for (i = 0; i < 4; i++)
    mynetwork[i] = (argc - 1) < i? 1: atoi(argv[i+1]);
  show_coor("mynet", 4, mynetwork, mynode);
  QMP_declare_logical_topology(mynetwork, 4);

  show_coor("allocated",
	    QMP_get_allocated_number_of_dimensions(),
	    QMP_get_allocated_dimensions(),
	    QMP_get_allocated_coordinates());
  show_coor("logical",
	    QMP_get_logical_number_of_dimensions(),
	    QMP_get_logical_dimensions(),
	    QMP_get_logical_coordinates());
  
  /* XXX */
  QMP_finalize_msg_passing();
  return 0;
}
