#include <mdwf.h>

static Up_project up_project_n[Q(DIM)] = {
  qx(proj_Ucg0plus),
  qx(proj_Ucg1plus),
  qx(proj_Ucg2plus),
  qx(proj_Ucg3plus)
};

static Down_project down_project_n[Q(DIM)] = {
  qx(proj_g0minus),
  qx(proj_g1minus),
  qx(proj_g2minus),
  qx(proj_g3minus)
};

static Up_project up_project_c[Q(DIM)] = {
  qx(proj_Ucg0minus),
  qx(proj_Ucg1minus),
  qx(proj_Ucg2minus),
  qx(proj_Ucg3minus)
};

static Down_project down_project_c[Q(DIM)] = {
  qx(proj_g0plus),
  qx(proj_g1plus),
  qx(proj_g2plus),
  qx(proj_g3plus)
};

/* BEGIN(prototypes) for mdwf.h */
/* END(prototypes) */	
	      
/*  r_x = a_x - Fxy b_y */
static void
qx(1_F)(struct Q(State) *state,
	struct eo_lattice *xy,
	struct Fermion *r_x,
	const struct Fermion *a_x,
	const struct SUn *U,
	const struct Fermion *b_y,
	long long *flops,
	long long *sent,
	long long *received)
{
    int Ls = state->Ls;
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    *flops += (up_project_n[i])(xy->send_up_buf[i],
					xy->send_up_size[i], Ls,
					xy->up_pack[i], U, b_y);
	if (xy->send_down_size[i])
	    *flops += (down_project_n[i])(xy->send_down_buf[i],
					  xy->send_down_size[i], Ls,
					  xy->down_pack[i], b_y);
    }
    if (xy->h_valid)
	QMP_start(xy->handle);

    *flops += qx(do_1mF)(r_x, 0, xy->body_size, Ls,
			 xy->body_neighbor, a_x, U, b_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    *flops += qx(do_1mF)(r_x, xy->body_size, xy->face_size, Ls,
			 xy->face_neighbor, a_x, U, b_y, xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;
}

/*  r_y = a_y - Fxy* b_x */
static void
qx(1_Fc)(struct Q(State) *state,
	 struct eo_lattice *yx,
	 struct Fermion *r_y,
	 const struct Fermion *a_y,
	 const struct SUn *U,
	 const struct Fermion *b_x,
	 long long *flops,
	 long long *sent,
	 long long *received)
{
    int Ls = state->Ls;
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (yx->send_up_size[i])
	    *flops += (up_project_c[i])(yx->send_up_buf[i],
					yx->send_up_size[i], Ls,
					yx->up_pack[i], U, b_x);
	if (yx->send_down_size[i])
	    *flops += (down_project_c[i])(yx->send_down_buf[i],
					  yx->send_down_size[i], Ls,
					  yx->down_pack[i], b_x);
    }
    if (yx->h_valid)
	QMP_start(yx->handle);

    *flops += qx(do_1mFc)(r_y, 0, yx->body_size, Ls,
			  yx->body_neighbor, a_y, U, b_x, NULL);

    if (yx->h_valid)
	QMP_wait(yx->handle);

    *flops += qx(do_1mFc)(r_y, yx->body_size, yx->face_size, Ls,
			  yx->face_neighbor, a_y, U, b_x, yx->receive_buf);
    *sent += yx->total_send;
    *received += yx->total_receive;
}

/*  r_y = 1/Ayy* Byy* Fxy* a_x */
static void
qx(1AcBcFc)(struct Q(State) *state,
	    struct eo_lattice *yx,
	    const struct Q(Parameters) *params,
	    struct Fermion *r_y,
	    const struct SUn *U,
	    const struct Fermion *a_x,
	    long long *flops,
	    long long *sent,
	    long long *received)
{
    int Ls = state->Ls;
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (yx->send_up_size[i])
	    *flops += (up_project_c[i])(yx->send_up_buf[i],
					yx->send_up_size[i], Ls,
					yx->up_pack[i], U, a_x);
	if (yx->send_down_size[i])
	    *flops += (down_project_c[i])(yx->send_down_buf[i],
					  yx->send_down_size[i], Ls,
					  yx->down_pack[i], a_x);
    }
    if (yx->h_valid)
	QMP_start(yx->handle);

    *flops += qx(do_1AcBcFc)(r_y, 0, yx->body_size, Ls,
			     params->AipTable, params->AimTable, params->BTable,
			     yx->body_neighbor, U, a_x, NULL);

    if (yx->h_valid)
	QMP_wait(yx->handle);

    *flops += qx(do_1AcBcFc)(r_y, yx->body_size, yx->face_size, Ls,
			     params->AipTable, params->AimTable, params->BTable,
			     yx->face_neighbor, U, a_x, yx->receive_buf);
    *sent += yx->total_send;
    *received += yx->total_receive;
}

/*  r_x = Bxx 1/Axx Fxy a_y */
static void
qx(B1AF)(struct Q(State) *state,
	 struct eo_lattice *xy,
	 const struct Q(Parameters) *params,
	 struct Fermion *r_x,
	 const struct SUn *U,
	 const struct Fermion *a_y,
	 long long *flops,
	 long long *sent,
	 long long *received)
{
    int Ls = state->Ls;
    int i;

    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    *flops += (up_project_n[i])(xy->send_up_buf[i],
					xy->send_up_size[i], Ls,
					xy->up_pack[i], U, a_y);
	if (xy->send_down_size[i])
	    *flops += (down_project_n[i])(xy->send_down_buf[i],
					  xy->send_down_size[i], Ls,
					  xy->down_pack[i], a_y);
    }
    if (xy->h_valid)
	QMP_start(xy->handle);

    *flops += qx(do_B1AF)(r_x, 0, xy->body_size, Ls,
			  params->BTable, params->AipTable, params->AimTable,
			  xy->body_neighbor, U, a_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    *flops += qx(do_B1AF)(r_x, xy->body_size, xy->face_size, Ls,
			  params->BTable, params->AipTable, params->AimTable,
			  xy->face_neighbor, U, a_y, xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;
}

/*  r_x = a_x - Bxx 1/Axx Fxy b_y */
static void
qx(1_B1AF)(struct Q(State) *state,
	   struct eo_lattice *xy,
	   const struct Q(Parameters) *params,
	   struct Fermion *r_x,
	   const struct Fermion *a_x,
	   const struct SUn *U,
	   const struct Fermion *b_y,
	   long long *flops,
	   long long *sent,
	   long long *received)
{
    int Ls = state->Ls;
    int i;
    
    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    *flops += (up_project_n[i])(xy->send_up_buf[i],
					xy->send_up_size[i], Ls,
					xy->up_pack[i], U, b_y);
	if (xy->send_down_size[i])
	    *flops += (down_project_n[i])(xy->send_down_buf[i],
					  xy->send_down_size[i], Ls,
					  xy->down_pack[i], b_y);
    }
    if (xy->h_valid)
	QMP_start(xy->handle);

    *flops += qx(do_1mB1AF)(r_x, 0, xy->body_size, Ls,
			    params->BTable, params->AipTable, params->AimTable,
			    xy->body_neighbor, a_x, U, b_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    *flops += qx(do_1mB1AF)(r_x, xy->body_size, xy->face_size, Ls,
			    params->BTable, params->AipTable, params->AimTable,
			    xy->face_neighbor, a_x, U, b_y, xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;
}

/*  r_x = a_x - Bxx 1/Axx Fxy b_y & |r_x|^2 */
static void
qx(1_B1AF_norm)(struct Q(State) *state,
		struct eo_lattice *xy,
		const struct Q(Parameters) *params,
		struct Fermion *r_x,
		double *norm,
		const struct Fermion *a_x,
		const struct SUn *U,
		const struct Fermion *b_y,
		long long *flops,
		long long *sent,
		long long *received)
{
    int Ls = state->Ls;
    int i;

    *norm = 0;
    
    for (i = 0; i < Q(DIM); i++) {
	if (xy->send_up_size[i])
	    *flops += (up_project_n[i])(xy->send_up_buf[i],
					xy->send_up_size[i], Ls,
					xy->up_pack[i], U, b_y);
	if (xy->send_down_size[i])
	    *flops += (down_project_n[i])(xy->send_down_buf[i],
					  xy->send_down_size[i], Ls,
					  xy->down_pack[i], b_y);
    }
    if (xy->h_valid)
	QMP_start(xy->handle);

    *flops += qx(do_1mB1AF_norm)(r_x, norm, 0, xy->body_size, Ls,
				 params->BTable,
				 params->AipTable, params->AimTable,
				 xy->body_neighbor, a_x, U, b_y, NULL);

    if (xy->h_valid)
	QMP_wait(xy->handle);

    *flops += qx(do_1mB1AF_norm)(r_x, norm, xy->body_size, xy->face_size, Ls,
				 params->BTable,
				 params->AipTable, params->AimTable,
				 xy->face_neighbor, a_x, U, b_y,
				 xy->receive_buf);
    *sent += xy->total_send;
    *received += xy->total_receive;
    QMP_sum_double(norm);
}


/* r_x = 1/Axx* Bxx* s_x */
static void
qx(1AcBc)(struct Q(State) *state,
	  const struct eo_lattice *yx,
	  const struct Q(Parameters) *params,
	  struct Fermion *r_x,
	  const struct Fermion *s_x,
	  long long *flops)
{
    *flops += qx(do_1AcBc)(r_x, yx->full_size, state->Ls,
			   params->AipTable, params->AimTable, params->BTable, 
			   s_x);
}

/* r_x = Bxx 1/Axx s_x */
static void
qx(B1A)(struct Q(State) *state,
	const struct eo_lattice *yx,
	const struct Q(Parameters) *params,
	struct Fermion *r_x,
	const struct Fermion *s_x,
	long long *flops)
{
    *flops += qx(do_B1A)(r_x, yx->full_size, state->Ls,
			 params->BTable, params->AipTable, params->AimTable,
			 s_x);
}

/* r_x = 1/Axx s_x */
static void
qx(1A)(struct Q(State) *state,
       const struct eo_lattice *yx,
       const struct Q(Parameters) *params,
       struct Fermion *r_x,
       const struct Fermion *s_x,
       long long *flops)
{
    *flops += qx(do_A_inverse)(r_x, yx->full_size, state->Ls,
			       params->AipTable, params->AimTable,
			       s_x);
}

/* r_x = 1/Bxx s_x */
static void
qx(1B)(struct Q(State) *state,
       const struct eo_lattice *yx,
       const struct Q(Parameters) *params,
       struct Fermion *r_x,
       const struct Fermion *s_x,
       long long *flops)
{
    *flops += qx(do_A_inverse)(r_x, yx->full_size, state->Ls,
			       params->BipTable, params->BimTable,
			       s_x);
}

/*   Mx = (1 - Foe* (1/Aoo* Boo* Feo* (1/Aee* Bee*)))
 *  r_x = Mx s_x
 */
static void
qx(Mx)(struct Q(State) *state,
       struct eo_lattice *xy,
       struct eo_lattice *yx,
       const struct Q(Parameters) *params,
       struct Fermion *r_x,
       const struct SUn *U,
       const struct Fermion *s_x,
       struct Fermion *t0_x,
       struct Fermion *t1_y,
       long long *flops,
       long long *sent,
       long long *received)
{
    qx(1AcBc)(state, yx, params, t0_x, s_x, flops);
    qx(1AcBcFc)(state, yx, params, t1_y, U, t0_x, flops, sent, received);
    qx(1_Fc)(state, xy, r_x, s_x, U, t1_y, flops, sent, received);
}

/*   M  = (1 - Bee 1/Aee Feo (Boo 1/Aoo Foe))
 * r_x = M s_x
 */
static void
qx(M)(struct Q(State) *state,
      struct eo_lattice *xy,
      struct eo_lattice *yx,
      const struct Q(Parameters) *params,
      struct Fermion *r_x,
      const struct SUn *U,
      const struct Fermion *s_x,
      struct Fermion *t0_y,
      long long *flops,
      long long *sent,
      long long *received)
{
    qx(B1AF)(state, yx, params, t0_y, U, s_x, flops, sent, received);
    qx(1_B1AF)(state, xy, params, r_x, s_x, U, t0_y, flops, sent, received);
}

static void
qx(M_norm)(struct Q(State) *state,
	   struct eo_lattice *xy,
	   struct eo_lattice *yx,
	   const struct Q(Parameters) *params,
	   struct Fermion *r_x,
	   double *norm,
	   const struct SUn *U,
	   const struct Fermion *s_x,
	   struct Fermion *t0_y,
	   long long *flops,
	   long long *sent,
	   long long *received)
{
    qx(B1AF)(state, yx, params, t0_y, U, s_x, flops, sent, received);
    qx(1_B1AF_norm)(state, xy, params, r_x, norm, s_x, U, t0_y,
		    flops, sent, received);
}

/* Solve MxMxx xi_x = varphi_x
 *
 *  x <- x0
 *  r <- b - M* M x, rho <- |r|^2
 *  p <- r
 *  k <- 0
 *  while (rho > epsilon) && (k < n) do
 *    z <- M p, zz <- |z|^2
 *    q <- M* z
 *    alpha <- rho/zz
 *    r <- r - alpha q, gamma <- |r|^2
 *    x <- x + alpha p
 *    if (gamma < epsilon) then
 *       rho <- gamma
 *       break
 *    end
 *    beta <- gamma / rho
 *    rho <- gamma
 *    p <- r + beta p
 *    k <- k + 1
 *  done
 *  return x, rho, k
 *
 */
static void
qx(cg_MxM)(struct Q(State) *state,
	   struct eo_lattice *xy,
	   struct eo_lattice *yx,
	   const struct Q(Parameters) *params,
	   struct Fermion *xi_x,
	   const struct SUn *U,
	   const struct Fermion *varphi_x,
	   const struct Fermion *x0_x,
	   struct Fermion *t1_y,
	   struct Fermion *t0_x,
	   struct Fermion *p_x,
	   struct Fermion *r_x,
	   struct Fermion *q_x,
	   struct Fermion *z_x,
	   int max_iter, int *out_iter,
	   double epsilon, double *out_eps,
	   long long *flops,
	   long long *sent,
	   long long *received)
{
    double alpha, beta, gamma, rho, zz;
    int Ls = state->Ls;
    int x_size = xy->full_size;
    int k;

    qx(f_copy)(x_size, Ls, xi_x, x0_x);
    qx(M)(state, xy, yx, params, q_x, U, x0_x, t1_y,
	  flops, sent, received);
    qx(Mx)(state, xy, yx, params, p_x, U, q_x, r_x, t1_y,
	   flops, sent, received);
    *flops += qx(f_add3)(x_size, Ls, r_x, varphi_x, -1.0, p_x);
    *flops += qx(f_norm)(x_size, Ls, &rho, r_x);
    /* XXX */ rho = 1.0;
    qx(f_copy)(x_size, Ls, p_x, r_x);
    
    for (k = 0; rho > epsilon && k < max_iter; k++) {
	qx(M_norm)(state, xy, yx, params, z_x, &zz, U, p_x, t1_y,
		   flops, sent, received);
	/* XXX */ zz = 1.0;
	qx(Mx)(state, xy, yx, params, q_x, U, z_x, t0_x, t1_y,
	       flops, sent, received);
	alpha = rho / zz;
	*flops += qx(f_add2_norm)(x_size, Ls, r_x, &gamma, -alpha, q_x);
	/* XXX */ gamma = 1.0;
	*flops += qx(f_add2)(x_size, Ls, xi_x, alpha, p_x);
	if (gamma < epsilon)
	    break;
	beta = gamma / rho;
	rho = gamma;
	*flops += qx(f_add2x)(x_size, Ls, p_x, beta, r_x);
    }
    *out_iter = k;
    *out_eps = rho;
}

int
QX(DDW_CG)(struct QX(Fermion) *result,
	   int *out_iterations,
	   double *out_epsilon,
	   const struct Q(Parameters) *params,
	   const struct QX(Fermion) *x_guess,
	   const struct QX(Gauge) *gauge,
	   const struct QX(Fermion) *rhs,
	   int max_iterations,
	   double epsilon)
{
    DECLARE_STATE;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    void *ptr;
    size_t ptr_size;
    void *temps;
    struct Fermion *t1_o;
    struct Fermion *t2_e;  /* also z_e */
    struct Fermion *phi_e; /* also q_e */
    struct Fermion *varphi_e;
    struct Fermion *xi_e;
    struct Fermion *r_e;
    struct Fermion *p_e;
    struct Fermion *t3_e;

    CHECK_ARG0(result);
    CHECK_ARGn(params, "DDW_CG");
    CHECK_ARGn(x_guess, "DDW_CG");
    CHECK_ARGn(gauge, "DDW_CG");
    CHECK_ARGn(rhs, "DDW_CG");

    if (q(setup_comm)(state, sizeof (REAL))) {
	return q(set_error)(state, 0, "DDWF_CG(): communication setup failed");
    }

    /* allocate locals */
    ptr = q(allocate_eo)(state, &ptr_size,
			 &temps, /* first local */
			 0, /* no header */
			 7, /* evens */
			 1, /* odds */
			 sizeof (REAL));
    if (ptr == 0) {
	return q(set_error)(state, 0, "DDWF_CG(): not enough memory");
    }
    /* setup other locals from temps */
    t1_o     = temps;
    t2_e     = temps = q(step_odd)(state, temps, sizeof (REAL));
    phi_e    = temps = q(step_even)(state, temps, sizeof (REAL));
    varphi_e = temps = q(step_even)(state, temps, sizeof (REAL));
    xi_e     = temps = q(step_even)(state, temps, sizeof (REAL));
    r_e      = temps = q(step_even)(state, temps, sizeof (REAL));
    p_e      = temps = q(step_even)(state, temps, sizeof (REAL));
    t3_e     = temps = q(step_even)(state, temps, sizeof (REAL));
    
    BEGIN_TIMING(state);
    /*  compute rhs for MxM */
    /*  t1_o = Boo 1/Aoo rhs_o */
    qx(B1A)(state, &state->odd, params, t1_o, rhs->odd, &flops);
    /*  t2_e = rhs_e - Feo t1_o */
    qx(1_F)(state, &state->even, t2_e, rhs->even, gauge->data, t1_o,
	    &flops, &sent, &received);
    /*  phi_e = Bee 1/Aee t2_e */
    qx(B1A)(state, &state->even, params, phi_e, t2_e, &flops);
    /*  varphi_e = Mx phi_e, t1_o, t2_e as temps */
    qx(Mx)(state, &state->even, &state->odd, params, varphi_e,
	   gauge->data, phi_e,
	   t2_e, t1_o, /* temps */
	   &flops, &sent, &received);

    /* solve MxM */
    /*   xi_e = 1/MxM varphi_e , with temps and results */
    qx(cg_MxM)(state, &state->even, &state->odd, params, xi_e,
	       gauge->data, varphi_e,
	       x_guess->even, /* XXX this is wrong */
	       t1_o, t3_e, p_e, r_e, phi_e, t2_e, /* temps */
	       max_iterations, out_iterations,
	       epsilon, out_epsilon,
	       &flops, &sent, &received);

    /* compute the other part of the result */
    /*  t1_o = rhs_o - Foe xi_e */
    qx(1_F)(state, &state->odd, t1_o, rhs->odd, gauge->data, xi_e,
	    &flops, &sent, &received);
    /*  result_o = 1/Aoo t1_o */
    qx(1A)(state, &state->odd, params, result->odd, t1_o, &flops);
    /*  result_e = 1/Bee xi_e */
    qx(1B)(state, &state->even, params, result->even, xi_e, &flops);

    END_TIMING(state, flops, sent, received);
    /* free locals */
    q(free)(state, ptr, ptr_size);
    return 0;
}
