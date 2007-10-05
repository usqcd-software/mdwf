#include <qmp.h>
#include <qop-mdwf3.h>

int self;

int
main(int argc, char *argv[])
{
  QMP_thread_level_t qt = QMP_THREAD_SINGLE;
  int net_dim;
  const int *network;
  const int *node;
  int i;
  int primary;

  if (QMP_init_msg_passing(&argc, &argv, qt, &qt) != QMP_SUCCESS) {
    fprintf(stderr, "QMP_init() failed\n");
    return 1;
  }

  self = QMP_get_node_number();
  primary = QMP_is_primary_node();
  net_dim = QMP_get_allocated_number_of_dimensions();
  network = QMP_get_allocated_dimensions();
  node = QMP_get_allocated_coordinates();
  printf("[%d] machine of dimension %d, primary? %d\n", self, net_dim, primary);
  printf("[%d] network:", self);
  for (i = 0; i < net_dim; i++)
    printf(" %d", network[i]);
  printf("\n");
  printf("[%d] node:", self);
  for (i = 0; i < net_dim; i++)
    printf(" %d", node[i]);
  printf("\n");

  
  /* XXX */
  QMP_finalize_msg_passing();
  return 0;
}
