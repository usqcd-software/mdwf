#ifndef MARK_0DAB983B_E5FA_4BCE_AE68_915CE2ABD086
#define MARK_0DAB983B_E5FA_4BCE_AE68_915CE2ABD086

void dot_fermion(double *v_r, double *v_i,
		 const struct QX(Fermion) *a,
		 const struct QX(Fermion) *b);
int op_A_even(struct Fermion *result,
	      const struct Q(Parameters) *params,
	      const struct Fermion *fermion);
int op_A_odd(struct Fermion *result,
	     const struct Q(Parameters) *params,
	     const struct Fermion *fermion);
int op_A(struct QX(Fermion) *result,
	 const struct Q(Parameters) *params,
	 const struct QX(Fermion) *fermion);
int op_Ax_even(struct Fermion *result,
	       const struct Q(Parameters) *params,
	       const struct Fermion *fermion);
int op_Ax_odd(struct Fermion *result,
	      const struct Q(Parameters) *params,
	      const struct Fermion *fermion);
int op_Ax(struct QX(Fermion) *result,
	  const struct Q(Parameters) *params,
	  const struct QX(Fermion) *fermion);
int op_A1_even(struct Fermion *result,
	       const struct Q(Parameters) *params,
	       const struct Fermion *fermion);
int op_A1_odd(struct Fermion *result,
	      const struct Q(Parameters) *params,
	      const struct Fermion *fermion);
int op_A1(struct QX(Fermion) *result,
	  const struct Q(Parameters) *params,
	  const struct QX(Fermion) *fermion);
int op_A1x_even(struct Fermion *result,
		const struct Q(Parameters) *params,
		const struct Fermion *fermion);
int op_A1x_odd(struct Fermion *result,
	       const struct Q(Parameters) *params,
	       const struct Fermion *fermion);
int op_A1x(struct QX(Fermion) *result,
	   const struct Q(Parameters) *params,
	   const struct QX(Fermion) *fermion);
int op_B_even(struct Fermion *result,
	      const struct Q(Parameters) *params,
	      const struct Fermion *fermion);
int op_B_odd(struct Fermion *result,
	     const struct Q(Parameters) *params,
	     const struct Fermion *fermion);
int op_B(struct QX(Fermion) *result,
	 const struct Q(Parameters) *params,
	 const struct QX(Fermion) *fermion);
int op_Bx_even(struct Fermion *result,
	       const struct Q(Parameters) *params,
	       const struct Fermion *fermion);
int op_Bx_odd(struct Fermion *result,
	      const struct Q(Parameters) *params,
	      const struct Fermion *fermion);
int op_Bx(struct QX(Fermion) *result,
	  const struct Q(Parameters) *params,
	  const struct QX(Fermion) *fermion);
int op_B1_even(struct Fermion *result,
	       const struct Q(Parameters) *params,
	       const struct Fermion *fermion);
int op_B1_odd(struct Fermion *result,
	      const struct Q(Parameters) *params,
	      const struct Fermion *fermion);
int op_B1(struct QX(Fermion) *result,
	  const struct Q(Parameters) *params,
	  const struct QX(Fermion) *fermion);
int op_B1x_even(struct Fermion *result,
		const struct Q(Parameters) *params,
		const struct Fermion *fermion);
int op_B1x_odd(struct Fermion *result,
	       const struct Q(Parameters) *params,
	       const struct Fermion *fermion);
int op_B1x(struct QX(Fermion) *result,
	   const struct Q(Parameters) *params,
	   const struct QX(Fermion) *fermion);
int op_F_even(struct Fermion *result_even,
	      struct Q(State) *state,
	      const struct SUn *U,
	      const struct Fermion *src_odd);
int op_F_odd(struct Fermion *result_odd,
	     struct Q(State) *state,
	     const struct SUn *U,
	     const struct Fermion *src_even);
int op_F(struct Q(Fermion) *result,
	 const struct Q(Gauge) *gauge,
	 const struct Q(Fermion) *source);
int op_Fx_even(struct Fermion *result_even,
	       struct Q(State) *state,
	       const struct SUn *U,
	       const struct Fermion *src_odd);
int op_Fx_odd(struct Fermion *result_odd,
	      struct Q(State) *state,
	      const struct SUn *U,
	      const struct Fermion *src_even);
int op_Fx(struct Q(Fermion) *result,
	  const struct Q(Gauge) *gauge,
	  const struct Q(Fermion) *source);
int op_A1xBx_even(struct Fermion *res,
		  const struct Q(Parameters) *params,
		  const struct Fermion *src);
int op_A1xBx_odd(struct Fermion *res,
		 const struct Q(Parameters) *params,
		 const struct Fermion *src);
int op_A1xBx(struct Q(Fermion) *res,
	     const struct Q(Parameters) *params,
	     const struct Q(Fermion) *src);
int op_ApF_odd(struct Fermion *r_odd,
	       struct Q(State) *state,
	       const struct Q(Parameters) *params,
	       const struct SUn *U,
	       const struct Fermion *a_odd,
	       const struct Fermion *a_even);
int op_1mFx_odd(struct Fermion *r_odd,
		struct Q(State) *state,
		const struct SUn *U,
		const struct Fermion *a_odd,
		const struct Fermion *a_even);
int op_1mBA1F_odd(struct Fermion *r_odd,
		  struct Q(State) *state,
		  const struct Q(Parameters) *params,
		  const struct SUn *U,
		  const struct Fermion *a_odd,
		  const struct Fermion *a_even);
int op_BA1Fn_odd(struct Fermion *r_odd,
		 double *norm,
		 struct Q(State) *state,
		 const struct Q(Parameters) *params,
		 const struct SUn *U,
		 const struct Fermion *a_odd,
		 const struct Fermion *a_even);
int op_A1xBxFx_odd(struct Fermion *r_odd,
		   struct Q(State) *state,
		   const struct Q(Parameters) *params,
		   const struct SUn *U,
		   const struct Fermion *a_even);
int op_BA1F_odd(struct Fermion *r_odd,
		struct Q(State) *state,
		const struct Q(Parameters) *params,
		const struct SUn *U,
		const struct Fermion *a_even);

extern Up_project up_project_n[Q(DIM)];
extern Up_project up_project_x[Q(DIM)];
extern Down_project down_project_n[Q(DIM)];
extern Down_project down_project_x[Q(DIM)];

#endif /* !defined(MARK_0DAB983B_E5FA_4BCE_AE68_915CE2ABD086) */
