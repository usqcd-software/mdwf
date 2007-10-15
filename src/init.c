#include <mdwf.h>

static void
set_send(struct eo_lattice *here,
	 const struct eo_lattice *there)
{
  int i;

  for (i = 0; i < Q(DIM); i++) {
    here->send_up_size[i] = there->receive_up_size[i];
    here->send_down_size[i] = there->receive_down_size[i];
   }
 }

struct p_layout {
  int *layout2vector;
  int *vector2layout;
};

struct layout {
  struct p_layout full;
  struct p_layout even;
  struct p_layout odd;
};

static int
alloc_p(struct p_layout *p,
	int full,
	int size,
	struct Q(State) *state)
{
  p->layout2vector = q(malloc)(state, size * sizeof (int));
  if (p->layout2vector == NULL)
    return 1;

  p->vector2layout = q(malloc)(state, full * sizeof (int));
  if (p->vector2layout == NULL)
    return 1;

  return 0;
}

static int
alloc_layout(struct layout *layout,
	     int full_size,
	     int even_size,
	     int odd_size,
	     struct Q(State) *state)
{
  if (alloc_p(&layout->full, full_size, full_size, state))
    return 1;
  if (alloc_p(&layout->even, full_size, even_size, state))
    return 1;
  if (alloc_p(&layout->odd, full_size, odd_size, state))
    return 1;
  return 0;
}

static int
v2l(const int lo[Q(DIM)], const int hi[Q(DIM)], const int x[Q(DIM)])
{
  int p, i;
  for (p = 0, i = 0; i < Q(DIM); i++) {
    p = p * (hi[i] - lo[i]) + x[i] - lo[i];
  }
  return p;
}

static void
l2v(int x[Q(DIM)], const int lo[Q(DIM)], const int hi[Q(DIM)], int p)
{
  int i;
  for (i = 0; i < Q(DIM); i++) {
    int d = hi[i] - lo[i];
    x[i] = p % d + lo[i];
    p = p / d;
  }
}

static void
walker(struct layout *layout,
       const int lo[Q(DIM)],
       int left[Q(DIM)],
       int right[Q(DIM)],
       const int hi[Q(DIM)],
       const int network[Q(DIM)],
       int *fi, int *ei, int *oi, int fp)
{
  int b, i, e;

  for (e = right[0] - left[0], b = 0, i = 1; i < Q(DIM); i++) {
    int g = right[i] - left[i];
    if (g > e) {
      e = g;
      b = i;
    }
  }
  if (e > 1) {
    int dl = left[b];
    int dr = right[b];
    right[b] = dl + e/2;
    walker(layout, lo, left, right, hi, network, fi, ei, oi, fp);
    left[b] = right[b];
    right[b] = dr;
    walker(layout, lo, left, right, hi, network, fi, ei, oi, fp);
    left[b] = dl;
  } else {
    int p, parity, face;

    for (parity = face = 0, i = 0; i < Q(DIM); i++) {
      parity += left[i];
      if (network[i] == 1)
	continue;
      if (left[i] == lo[i] || left[i] + 1 == hi[i])
	face = 1;
    }
    if (face != fp)
      return;

    p = v2l(lo, hi, left);
    layout->full.vector2layout[p] = *fi;
    layout->full.layout2vector[*fi++] = p;
    if (parity & 1) {
      layout->odd.vector2layout[p] = *oi;
      layout->odd.layout2vector[*oi++] = p;      
    } else {
      layout->even.vector2layout[p] = *ei;
      layout->even.layout2vector[*ei++] = p;
    }
  }
}
       

static void
build_layout(struct sublattice *sublattice,
	     const int network[Q(DIM)],
	     struct layout *layout)
	     
{
  int f = 0;
  int e = 0;
  int o = 0;
  int i;
  const int *lo = sublattice->lo;
  const int *hi = sublattice->hi;
  int left[Q(DIM)];
  int right[Q(DIM)];

  for (i = 0; i < Q(DIM); i++) {
    left[i] = lo[i];
    right[i] = hi[i];
  }

  walker(layout, lo, left, right, hi, network, &f, &e, &o, 0);
  walker(layout, lo, left, right, hi, network, &f, &e, &o, 1);
}

static void
free_p(struct p_layout *p, int full, int size, struct Q(State) *state)
{
  if (p->layout2vector)
    q(free)(state, p->layout2vector, size * sizeof (int));
  p->layout2vector = NULL;

  if (p->vector2layout)
    q(free)(state, p->vector2layout, full * sizeof (int));
  p->vector2layout= NULL;
}

static void
free_layout(struct layout *layout,
	    int full_size,
	    int even_size,
	    int odd_size,
	    struct Q(State) *state)
{
  free_p(&layout->full, full_size, full_size, state);
  free_p(&layout->even, full_size, even_size, state);
  free_p(&layout->odd, full_size, odd_size, state);
}

static int
eo_alloc(struct Q(State) *state,
	  struct eo_lattice *eo)
{
  int i;
  int up_pack_size;
  int neighbor_size;
  const char *error = NULL;
  
  q(sizeof_up_pack)(&up_pack_size);
  q(sizeof_neighbor)(&neighbor_size);
  for (i = 0; i < Q(DIM); i++) {
    if (eo->send_up_size[i]) {
      eo->up_pack[i] = q(malloc)(state, eo->send_up_size[i] * up_pack_size);
      if (eo->up_pack[i] == NULL) {
	error = "Not enough space for E/O up_pack";
	goto failed;
      }
    }
    if (eo->send_down_size[i]) {
      eo->down_pack[i] = q(malloc)(state, eo->send_down_size[i] * sizeof (int));
      if (eo->down_pack[i] == NULL) {
	error = "Now enough space for E/O down/pack";
	goto failed;
      }
    }
  }
  if (eo->full_size) {
    eo->body_neighbor = q(malloc)(state, eo->full_size * neighbor_size);
    if (eo->body_neighbor == NULL) {
      error = "Not enough space for E/O body_neighbor";
      goto failed;
    }
    eo->face_neighbor = (struct neighbor *)((char *)eo->body_neighbor
					    + eo->body_size * neighbor_size);
  }
  return 0;

 failed:
  for (i = 0; i < Q(DIM); i++) {
    if (eo->up_pack[i])
      q(free)(state, eo->up_pack[i], eo->send_up_size[i] * up_pack_size);
    eo->up_pack[i] = NULL;
    if (eo->down_pack[i])
      q(free)(state, eo->down_pack[i], eo->send_down_size[i] * sizeof (int));
    eo->down_pack[i] = NULL;
  }
  if (eo->body_neighbor)
    q(free)(state, eo->body_neighbor, eo->body_size * neighbor_size);
  eo->body_neighbor = 0;
  eo->face_neighbor = 0;
  
  q(set_error)(state, 1, error);
  return 1;
}

static void
build_local_neighbors(struct eo_lattice *eo,
		      struct Q(State) *state,
		      const int network[Q(DIM)],
		      const int *full_v2l,
		      const int *oe_v2l)
{
  int p, d;
  int x[Q(DIM)];
  int la, lb;
  int full_size = eo->full_size;
  int mask;
  int u_up;
  int u_down[Q(DIM)];
  int f_up[Q(DIM)];
  int f_down[Q(DIM)];
  const int *lo = state->sublattice.lo;
  const int *hi = state->sublattice.hi;
  
  for (p = 0; p < full_size; p++) {
    la = eo->layout2vector[p];
    l2v(x, lo, hi, la);
    mask = 0;
    u_up = full_v2l[la] * Q(DIM);
    for (d = 0; d < Q(DIM); d++) {
      x[d] += 1;
      if ((network[d] > 1) && (x[d] == hi[d])) {
	mask |= 1 << d;
	f_up[d] = -1;
      } else {
	if (x[d] == hi[d]) x[d] = 0;
	f_up[d] = oe_v2l[v2l(lo, hi, x)];
	if (x[d] == 0) x[d] = hi[d];
      }
      x[d] -= 2;
      if ((network[d] > 1) && (x[d] < state->sublattice.lo[d])) {
	mask |= 1 << (d + Q(DIM));
	f_down[d] = -1;
	u_down[d] = -1;
      } else {
	if (x[d] < 0) x[d] = hi[d] - 1;
	lb = v2l(lo, hi, x);
	f_down[d] = oe_v2l[lb];
	u_down[d] = full_v2l[lb] * Q(DIM) + d;
	if (x[d] == hi[d] - 1) x[d] = -1;
      }
      x[d] += 1;
    }
    q(set_neighbor)(eo->body_neighbor, p, mask, f_up, u_up, f_down, u_down);
  }
}

static int
eo_init(struct Q(State) *state,
	struct eo_lattice *this,
	struct eo_lattice *other)
{
  /* XXX */
  return 0;
}

int
Q(init)(struct Q(State) **state_ptr,
	const int lattice[Q(DIM)+1],
	const int network[Q(DIM)],
	const int node[Q(DIM)],
	int master_p,
	void (*sublattice)(int lo[],
			   int hi[],
			   const int node[],
			   void *env),
	void *env)
{
  struct Q(State) *state;
  int p, q, i;
  int dx[Q(DIM)], x[Q(DIM)];
  int volume;
  struct layout layout;
  const char *error;
#define INIT_ERROR(msg) do { error = msg; goto failed; } while (0)

  if (state_ptr == NULL)
    return 1;

  *state_ptr = NULL;
  state = q(malloc)(NULL, sizeof (struct Q(State)));
  if (state == NULL)
    return 1;

  memset(&layout, 0, sizeof (struct layout));
  memset(state, 0, sizeof (struct Q(State)));
  *state_ptr = state;
  
  state->used = 1;
  state->Ls = lattice[Q(DIM)];
  state->allocated = state->max_allocated = sizeof (struct Q(State));
  sublattice(state->sublattice.lo, state->sublattice.hi, node, env);
  for (volume = 1, i = 0; i < Q(DIM); i++) {
    dx[i] = state->sublattice.hi[i] - state->sublattice.lo[i];
    volume *= dx[i];
    if (lattice[i] % 2)
      INIT_ERROR("Odd lattice size");
    state->lattice[i] = lattice[i];
    if (network[i] < 1)
      INIT_ERROR("Zero size network");
    state->network[i] = network[i];
    if (node[i] < 0)
      INIT_ERROR("Bad network coordinates");
    state->node[i] = node[i];
  }
  state->volume = volume;
  state->master_p = master_p;
  for (p = 0; p < volume; p++) {
    struct eo_lattice *here;
    int parity = 0;
    int face = 0;
    int face_up[Q(DIM)];
    int face_down[Q(DIM)];

    for (i = 0; i < Q(DIM); i++) {
      face_up[i] = face_down[i] = 0;
    }
    for (q = p, i = 0; i < Q(DIM); i++) {
      x[i] = q % dx[i] + state->sublattice.lo[i];
      q = q / dx[i];
      parity += x[i];
      if (x[i] == state->sublattice.lo[i]) {
	face_down[i] = 1;
	face = 1;
      }
      if (x[i] + 1 == state->sublattice.hi[i]) {
	face_up[i] = 1;
	face = 1;
      }
    }
    here = parity == 0? &state->even: &state->odd;
    here->full_size++;
    if (face)
      here->face_size++;
    else
      here->body_size++;
    for (i = 0; i < Q(DIM); i++) {
      here->receive_down_size[i] += face_down[i];
      here->receive_up_size[i] += face_down[i];
    }
  }
  set_send(&state->even, &state->odd);
  set_send(&state->odd, &state->even);

  if (eo_alloc(state, &state->even))
    INIT_ERROR(NULL);
  if (eo_alloc(state, &state->odd))
    INIT_ERROR(NULL);

  if (alloc_layout(&layout,
		   state->volume,
		   state->even.full_size,
		   state->odd.full_size,
		   state))
    INIT_ERROR("Not enough memory for layout");

  build_layout(&state->sublattice, state->network, &layout);

  state->even.layout2vector = layout.even.layout2vector;
  layout.even.layout2vector = NULL;
  state->odd.layout2vector = layout.odd.layout2vector;
  layout.odd.layout2vector = NULL;
  state->layout2vector = layout.full.layout2vector;
  layout.full.layout2vector = NULL;

  build_local_neighbors(&state->even, state, network,
			layout.full.vector2layout,
			layout.odd.vector2layout);
  build_local_neighbors(&state->odd, state, network,
			layout.full.vector2layout,
			layout.even.vector2layout);
  /* XXX */

  free_layout(&layout,
	      state->volume,
	      state->even.full_size,
	      state->odd.full_size,
	      state);

  /* XXX */

  if (eo_init(state, &state->even, &state->odd))
    INIT_ERROR(NULL);
  if (eo_init(state, &state->odd, &state->even))
    INIT_ERROR(NULL);

  /* XXX */

  return 0;

 failed:
  free_layout(&layout,
	      state->volume,
	      state->even.full_size,
	      state->odd.full_size,
	      state);
  if (error == NULL)
    error = state->error;
  q(cleanup_state)(state);
  q(set_error)(state, 1, error);

  return 1;
}
