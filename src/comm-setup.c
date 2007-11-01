#include <mdwf.h>

static int
eo_comm(struct eo_lattice *eo, struct Q(State) *state, int real_size)
{
  int d, k, ks, kr;
  int pfL = eo->Ls * Q(COLORS) * Q(PROJECTED_FERMION_DIM) * 2 * real_size;
  int size;

  if (eo->real_size == real_size)
    return 0;

  xprint("eo_comm(%p, %p, %d), receive_buf=%p",
	 eo, state, real_size, eo->receive_buf);

  /* This version bundles all communication together. Seems to work, but may
   * be suboptimal.
   */
  
  eo->total_send = 0;
  eo->total_receive = 0;
  eo->real_size = real_size;
  for (d = 0, k = ks = kr = 0; d < Q(DIM); d++) {
    if (eo->send_up_size[d]) {
      size = eo->send_up_size[d] * pfL;
      eo->mem[k] = QMP_allocate_aligned_memory(size,
					       CACHE_LINE_SIZE,
					       QMP_MEM_DEFAULT);
      if (eo->mem[k] == NULL)
	return 1;
      eo->mem_count = k + 1;
      eo->send_up_buf[d] = QMP_get_memory_pointer(eo->mem[k]);
      eo->mh[k] = QMP_declare_msgmem(eo->send_up_buf[d], size);
      if (eo->mh[k] == NULL)
	return 1;
      eo->mh_count = k + 1;
      eo->hs_up[d] = QMP_declare_send_relative(eo->mh[k], d, +1, 0);
      if (eo->hs_up[d] == NULL)
	return 1;
      eo->hsv[ks] = eo->hs_up[d];
      eo->hsv_count = ks + 1;
      k++;
      ks++;
      eo->total_send += size;
    }
    if (eo->send_down_size[d]) {
      size = eo->send_down_size[d] * pfL;
      eo->mem[k] = QMP_allocate_aligned_memory(size,
					       CACHE_LINE_SIZE,
					       QMP_MEM_DEFAULT);
      if (eo->mem[k] == NULL)
	return 1;
      eo->mem_count = k + 1;
      eo->send_down_buf[d] = QMP_get_memory_pointer(eo->mem[k]);
      eo->mh[k] = QMP_declare_msgmem(eo->send_down_buf[d], size);
      if (eo->mh[k] == NULL)
	return 1;
      eo->mh_count = k + 1;
      eo->hs_down[d] = QMP_declare_send_relative(eo->mh[k], d, -1, 0);
      if (eo->hs_down[d] == NULL)
	return 1;
      eo->hsv[ks] = eo->hs_down[d];
      eo->hsv_count = ks + 1;
      k++;
      ks++;
      eo->total_send += size;
    }
    if (eo->receive_up_size[d]) {
      size = eo->receive_up_size[d] * pfL;
      eo->mem[k] = QMP_allocate_aligned_memory(size,
					       CACHE_LINE_SIZE,
					       QMP_MEM_DEFAULT);
      if (eo->mem[k] == NULL)
	return 1;
      eo->mem_count = k + 1;
      eo->receive_buf[d] = QMP_get_memory_pointer(eo->mem[k]);
      eo->mh[k] = QMP_declare_msgmem(eo->receive_buf[d], size);
      if (eo->mh[k] == NULL)
	return 1;
      eo->mh_count = k + 1;
      eo->hrt[kr] = QMP_declare_receive_relative(eo->mh[k], d, +1, 0);
      if (eo->hrt[kr] == NULL)
	return 1;
      eo->hrt_count = kr + 1;
      k++;
      kr++;
      eo->total_receive += size;
    }
    if (eo->receive_down_size[d]) {
      size = eo->receive_down_size[d] * pfL;
      eo->mem[k] = QMP_allocate_aligned_memory(size,
					       CACHE_LINE_SIZE,
					       QMP_MEM_DEFAULT);
      if (eo->mem[k] == NULL)
	return 1;
      eo->mem_count = k + 1;
      eo->receive_buf[d + Q(DIM)] = QMP_get_memory_pointer(eo->mem[k]);
      eo->mh[k] = QMP_declare_msgmem(eo->receive_buf[d + Q(DIM)], size);
      if (eo->mh[k] == NULL)
	return 1;
      eo->mh_count = k + 1;
      eo->hrt[kr] = QMP_declare_receive_relative(eo->mh[k], d, -1, 0);
      if (eo->hrt[kr] == NULL)
	return 1;
      eo->hrt_count = kr + 1;
      k++;
      kr++;
      eo->total_receive += size;
    }
  }

  xprint("comm_setup(eo=%p, receive_buf[]={%p %p %p %p %p %p %p %p}",
	 eo,
	 eo->receive_buf[0],
	 eo->receive_buf[1],
	 eo->receive_buf[2],
	 eo->receive_buf[3],
	 eo->receive_buf[4],
	 eo->receive_buf[5],
	 eo->receive_buf[6],
	 eo->receive_buf[7]);

  return 0;
}

int
q(setup_comm)(struct Q(State) *state, int real_size)
{
  xprint("setup_comm(%d)", real_size);
  
  if (state->real_size == real_size)
    return 0;

  q(free_comm)(state);
  state->real_size = real_size;
  if (eo_comm(&state->even, state, real_size))
    return q(free_comm)(state);
  if (eo_comm(&state->odd, state, real_size))
    return q(free_comm)(state);

  xprint("setup_comm(%d) done", real_size);

  return 0;
}
