#include <mdwf.h>

void
q(l2v)(int x[Q(DIM)], const struct local *local, int la)
{
  int d;
  for (d = Q(DIM); d--;) {
    x[d] = la % local->dx[d] + local->lo[d];
    la = la / local->dx[d];
  }
}

int
q(v2l)(const int x[Q(DIM)], const struct local *local)
{
  int p, d;

  for (p = 0, d = 0; d < Q(DIM); d++) {
    p = p * local->dx[d] + x[d] - local->lo[d];
  }
  return p;
}

static void
build_packs(struct eo_lattice   *eo,
	    struct eo_lattice   *oe,
	    struct Q(State)     *state)
{
  int            d;
  int            z_up[Q(DIM)];
  int            z_down[Q(DIM)];
  int            p;
  int            x_bsize = oe->body_size;
  int            x_fsize = oe->full_size;
  struct local  *local = eo->local;
  int            x[Q(DIM)];
  
  for (d = 0; d < Q(DIM); d++) {
    z_up[d] = z_down[d] = 0;
  }
  for (p = x_bsize; p < x_fsize; p++) {
    int la = oe->lx2v[p];
    q(l2v)(x, oe->local, la);
    for (d = 0; d < Q(DIM); d++) {
      if (state->network[d] == 1)
	continue;
      if (x[d] == local->lo[d]) {
	q(put_down_pack)(eo->down_pack[d], z_down[d]++, p);
      }
      if (x[d] + 1 == local->hi[d]) {
	q(put_up_pack)(eo->up_pack[d], z_up[d]++,
		       p, state->v2lx[la] * Q(DIM) + d);
      }
    }
  }
}

static void
build_local_neighbors(struct eo_lattice        *eo,
		      const struct eo_lattice  *oe,
		      const struct Q(State)    *state)
{
  int                 p;
  int                 d;
  int                 x[Q(DIM)];
  int                 la;
  const struct local *local = eo->local;
  int                 fsize = eo->full_size;
  int                 mask;
  int                 f_up[Q(DIM)];
  int                 u_up;
  int                 f_down[Q(DIM)];
  int                 u_down[Q(DIM)];
  
  for (p = 0; p < fsize; p++) {
    la = eo->lx2v[p];
    q(l2v)(x, local, la);
    mask = 0;
    u_up = state->v2lx[la] * Q(DIM);
    for (d = 0; d < Q(DIM); d++) {
      x[d] += 1;
      if ((state->network[d] > 1) && (x[d] == local->hi[d])) {
	mask |= 1 << d;
	f_up[d] = -1;
      } else {
	if (x[d] == local->hi[d]) x[d] = 0;
	la = q(v2l)(x, local);
	f_up[d] = oe->v2lx[la];
	if (x[d] == 0) x[d] = local->hi[d];
      }
      x[d] -= 2;
      if ((state->network[d] > 1) && x[d] < local->lo[d]) {
	mask |= 1 << (d + Q(DIM));
	f_down[d] = -1;
	u_down[d] = -1;
      } else {
	if (x[d] < 0) x[d] = local->hi[d] - 1;
	la = q(v2l)(x, local);
	f_down[d] = oe->v2lx[la];
	u_down[d] = state->v2lx[la] * Q(DIM) + d;
	if (x[d] == local->hi[d] - 1) x[d] = -1;
      }
      x[d] += 1;
    }
    q(put_neighbor)(eo->body_neighbor, p, mask, f_up, u_up, f_down, u_down);
  }
}

static void
walker(struct Q(State)  *state,
       int               left[Q(DIM)],
       int               right[Q(DIM)],
       int               face_p,
       int              *index,
       int              *e_idx,
       int              *o_idx)
{
  int e, b, d;

  for (e = right[0] - left[0], b = 0, d = 1; d < Q(DIM); d++) {
    int g = right[d] - left[d];
    if (g < e)
      continue;
    e = g;
    b = d;
  }
  if (e > 1) {
    int dl = left[b];
    int dr = right[b];
    right[b] = dl + e/2;
    walker(state, left, right, face_p, index, e_idx, o_idx);
    left[b] = right[b];
    right[b] = dr;
    walker(state, left, right, face_p, index, e_idx, o_idx);
    left[b] = dl;
  } else {
    int la = q(v2l)(left, &state->local);
    int face = 0;
    int parity = 0;

    for (d = 0; d < Q(DIM); d++) {
      parity += left[d];
      if (state->network[d] > 1) {
	if ((left[d] == state->local.lo[d]) ||
	    (left[d] + 1 == state->local.hi[d]))
	  face = 1;
      }
    }
    if (face != face_p)
      return;

    state->v2lx[la] = *index;
    state->lx2v[*index] = la;
    *index += 1;
    if (parity & 1) {
      state->odd.v2lx[la] = *o_idx;
      state->odd.lx2v[*o_idx] = la;
      *o_idx += 1;
    } else {
      state->even.v2lx[la] = *e_idx;
      state->even.lx2v[*e_idx] = la;
      *e_idx += 1;
    }
  }
}

static void
build_layout(struct Q(State) *state)
{
  
  int  d;
  int  left[Q(DIM)];
  int  right[Q(DIM)];
  int  idx = 0;
  int  e_idx = 0;
  int  o_idx = 0;
  

  for (d = 0; d < Q(DIM); d++) {
    left[d] = state->local.lo[d];
    right[d] = state->local.hi[d];
  }
  walker(state, left, right, 0, &idx, &e_idx, &o_idx);
  walker(state, left, right, 1, &idx, &e_idx, &o_idx);
}

static char *
eo_init(struct eo_lattice        *eo,
	const struct eo_lattice  *oe,
	struct Q(State)          *state)
{
  int ns = q(sizeof_neighbor)(eo->full_size);
  int bs = q(sizeof_neighbor)(eo->body_size);
  int d;

#define CHECK(cond, msg) do { if (cond) return msg; } while (0)

  eo->Ls = state->Ls;
  eo->local = &state->local;
  eo->lx2v = q(malloc)(state, eo->full_size * sizeof(int));
  CHECK(eo->lx2v == NULL, "E/O->lx2v allocation failed");
  eo->v2lx = q(malloc)(state, state->volume * sizeof(int));
  CHECK(eo->v2lx == NULL, "E/O->v2lx allocation failed");
  eo->body_neighbor = q(malloc) (state, ns);
  CHECK(eo->body_neighbor == NULL, "E/O neighbor allocation failed");
  eo->face_neighbor = (struct neighbor*)(((char *)eo->body_neighbor) + bs);
  for (d = 0; d < Q(DIM); d++) {
    eo->send_up_size[d] = oe->receive_up_size[d];
    if (eo->send_up_size[d]) {
      int us = q(sizeof_up_pack)(eo->send_up_size[d]);
      eo->up_pack[d] = q(malloc)(state, us);
      CHECK(eo->up_pack[d] == NULL, "Not enough space for eo.up_pack[i]");
    }
    eo->send_down_size[d] = oe->receive_down_size[d];
    if (eo->send_down_size[d]) {
      int ds = q(sizeof_down_pack)(eo->send_down_size[d]);
      eo->down_pack[d] = q(malloc)(state, ds);
      CHECK(eo->down_pack[d] == NULL, "Not enough space for eo.down_pack[i]");
    }
  }
  return 0;
#undef CHECK
}

static char *
state_init(struct Q(State)  *state,
	   const int         lattice[Q(DIM) + 1],
	   const int         network[Q(DIM)],
	   const int         node[Q(DIM)],
	   int               master_p,
	   void            (*local)(int lo[],
				    int hi[],
				    const int node[],
				    void *env),
	   void             *env)
{
  char *msg;
  int v, d, la, parity, face;
  struct eo_lattice *eo;
  int x[Q(DIM)];
#define CHECK(cond,message) do { if (cond) return message; } while (0)

  memset(state, 0, sizeof (struct Q(State)));
  state->Ls = lattice[Q(DIM)];
  state->master_p = master_p;
  local(state->local.lo, state->local.hi, node, env);
  for (v = 1, d = 0; d < Q(DIM); d++) {
    state->network[d] = network[d];
    state->node[d] = node[d];
    state->lattice[d] = lattice[d];
    state->local.dx[d] = state->local.hi[d] - state->local.lo[d];
    v *= state->local.dx[d];
  }
  state->volume = v;

  for (la = 0; la < v; la++) {
    q(l2v)(x, &state->local, la);
    for (parity = 0, d = 0; d < Q(DIM); d++) {
      parity += x[d];
    }
    eo = (parity & 1) ? &state->odd : &state->even;
    for (face = 0, d = 0; d < Q(DIM); d++) {
      if (network[d] > 1) {
	if (x[d] == state->local.lo[d]) {
	  face = 1;
	  eo->receive_down_size[d]++;
	}
	if (x[d] == state->local.hi[d] - 1) {
	  face = 1;
	  eo->receive_up_size[d]++;
	}
      }
    }
    if (face)
      eo->face_size++;
    else
      eo->body_size++;
    eo->full_size++;
  }
  
  msg = eo_init(&state->even, &state->odd, state);
  CHECK(msg, msg);
  msg = eo_init(&state->odd, &state->even, state);
  CHECK(msg, msg);

  state->v2lx = q(malloc) (state, state->volume * sizeof(int));
  CHECK(state->v2lx == 0, "Not enough memory for state->vector2layout");

  state->lx2v = q(malloc) (state, state->volume * sizeof(int));
  CHECK(state->lx2v == 0, "Not enough memory for state->layout2vector");

  build_layout(state);

  build_local_neighbors(&state->even, &state->odd, state);
  build_local_neighbors(&state->odd, &state->even, state);
  build_packs(&state->even, &state->odd, state);
  build_packs(&state->odd, &state->even, state);

  return 0;
#undef CHECK
}

static void
eo_patch_up(struct eo_lattice *eo,
	    const struct Q(State) *state,
	    const struct eo_lattice *x_eo,
	    const struct eo_lattice *x_oe,
	    const struct Q(State) *x_state,
	    const int lattice[Q(DIM)+1],
	    int dim)
{
  int p, b, la, dp;
  int down_size = x_eo->send_down_size[dim];
  const struct local *local = &state->local;
  const struct local *x_local = &x_state->local;
  int x[Q(DIM)];

  for (p = 0; p < down_size; p++) {
    dp = q(get_down_pack_f)(x_eo->down_pack[dim], p);
    q(l2v)(x, x_local, x_oe->lx2v[dp]);
    x[dim]--;
    if (x[dim] < 0) x[dim] = lattice[dim] - 1;
    la = q(v2l)(x, local);
    b = eo->v2lx[la];
    q(fix_neighbor_f_up)(eo->body_neighbor, b, p, dim);
  }
}

static void
eo_patch_down(struct eo_lattice *eo,
	      const struct Q(State) *state,
	      const struct eo_lattice *x_eo,
	      const struct eo_lattice *x_oe,
	      const struct Q(State) *x_state,
	      const int lattice[Q(DIM)+1],
	      int dim)
{
  int p, b, la, up;
  int up_size = x_eo->send_up_size[dim];
  const struct local *local = &state->local;
  const struct local *x_local = &x_state->local;
  int x[Q(DIM)];

  for (p = 0; p < up_size; p++) {
    up = q(get_up_pack_f)(x_eo->up_pack[dim], p);
    q(l2v)(x, x_local, x_oe->lx2v[up]);
    x[dim]++;
    if (x[dim] == lattice[dim]) x[dim] = 0;
    la = q(v2l)(x, local);
    b = eo->v2lx[la];
    q(fix_neighbor_f_down)(eo->body_neighbor, b, p, dim);
  }
}

static char *
patch_boundary(struct Q(State) *state,
	       const int        lattice[Q(DIM) + 1],
	       const int        network[Q(DIM)],
	       const int        node[Q(DIM)],
	       int              master_p,
	       void           (*local) (int lo[], int hi[], const int node[],
					void *env),
	       void             *env)
{
  int d;
  int n[Q(DIM)];
  struct Q(State) x_state;
  char *status;
  
#define CHECK(cond,msg) do { if (cond) return msg;} while (0)
  for (d = 0; d < Q(DIM); d++)
    n[d] = node[d];

  for (d = 0; d < Q(DIM); d++) {
    if (network[d] == 1)
      continue;
    
    n[d]--;
    if (n[d] < 0) n[d] = network[d] - 1;

    status = state_init(&x_state, lattice, network, n, master_p, local, env);
    CHECK(status, "neighbor state init failed");
    eo_patch_down(&state->even, state,
		  &x_state.even, &x_state.odd, &x_state, lattice, d);
    eo_patch_down(&state->odd, state,
		  &x_state.odd, &x_state.even, &x_state, lattice, d);
    q(cleanup_state)(&x_state);

    n[d]++;
    if (n[d] == network[d]) n[d] = 0;

    n[d]++;
    if (n[d] == network[d]) n[d] = 0;

    status = state_init(&x_state, lattice, network, n, master_p, local, env);
    CHECK(status, "neighbor state init failed");
    eo_patch_up(&state->even, state,
		&x_state.even, &x_state.odd, &x_state, lattice, d);
    eo_patch_up(&state->odd, state,
		&x_state.odd, &x_state.even, &x_state, lattice, d);
    q(cleanup_state)(&x_state);

    n[d]--;
    if (n[d] < 0) n[d] = network[d] - 1;
  }
#undef CHECK
  return NULL;
}

int
Q(init)(struct Q(State) **state_ptr,
	const int lattice[Q(DIM)+1],
	const int network[Q(DIM)],
	const int node[Q(DIM)],
	int master_p,
	void (*local)(int lo[],
		      int hi[],
		      const int node[],
		      void *env),
	void *env)
{
  char *status;
  struct Q(State) *state;
#define CHECK(cond,message) do  {\
       if (cond) { status = message; goto error; } } while (0)

  if (state_ptr == 0)
    return 1;
  
  *state_ptr = state = q(malloc) (NULL, sizeof(struct Q(State)));
  if (state == 0)
    return 1;
  
  CHECK(lattice[0] % 1 != 0, "Lattice dimension X is not even");
  CHECK(lattice[1] % 1 != 0, "Lattice dimension Y is not even");
  CHECK(lattice[2] % 1 != 0, "Lattice dimension Z is not even");
  CHECK(lattice[3] % 1 != 0, "Lattice dimension T is not even");
  CHECK(lattice[0] < network[0], "Network is too large in X");
  CHECK(lattice[1] < network[1], "Network is too large in Y");
  CHECK(lattice[2] < network[2], "Network is too large in Z");
  CHECK(lattice[3] < network[3], "Network is too large in T");
  CHECK(network[0] < 1, "Network is too small in X");
  CHECK(network[1] < 1, "Network is too small in Y");
  CHECK(network[2] < 1, "Network is too small in Z");
  CHECK(network[3] < 1, "Network is too small in T");
  CHECK(node[0] >= network[0], "Node address too large in X");
  CHECK(node[1] >= network[1], "Node address too large in Y");
  CHECK(node[2] >= network[2], "Node address too large in Z");
  CHECK(node[3] >= network[3], "Node address too large in T");
  CHECK(node[0] < 0, "Node address too small in X");
  CHECK(node[1] < 0, "Node address too small in Y");
  CHECK(node[2] < 0, "Node address too small in Z");
  CHECK(node[3] < 0, "Node address too small in T");
  
  status = state_init(state, lattice, network, node, master_p, local, env);
  CHECK(status, status);
  
  status = patch_boundary(state, lattice, network, node, master_p, local, env);
  CHECK(status, status);
  
  q(free)(state, state->v2lx, state->volume * sizeof(int));
  state->v2lx = NULL;
  q(free)(state, state->even.v2lx, state->volume * sizeof (int));
  state->even.v2lx = NULL;
  q(free)(state, state->odd.v2lx, state->volume * sizeof (int));
  state->odd.v2lx = NULL;

  return 0;

error:
  q(cleanup_state)(state);
  q(set_error)(state, 1, status);
  return 1;
#undef CHECK
}
