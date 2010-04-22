#include <mdwf.h>

/* NOTE: This is probably true only for Shamir case */

int
QX(axial_current)(void                      (*writer)(const int pos[Q(DIM)],
                                                      int dir,
                                                      double val_re,
                                                      double val_im,
                                                      void *env),
                  void                       *env,
                  const struct QX(Fermion)   *fermion,
                  const struct QX(Gauge)     *gauge)
{
    DECLARE_STATE;
    long long flops = 0;
    long long sent = 0;
    long long received = 0;
    
    CHECK_ARG0(fermion);
    CHECK_ARGn(gauge, "axial_current");
    
    if (q(setup_comm)(state, sizeof (REAL)))
      return q(set_error)(state, 0,
                          "axial_current(): communication setup failed");
    
    BEGIN_TIMING(state);
    qx(op_axial_current)(writer, env,
                         &state->even, &state->odd, gauge->data,
                         fermion->even, fermion->odd,
                         &flops, &sent, &received, state->node);
    qx(op_axial_current)(writer, env,
                         &state->odd, &state->even, gauge->data,
                         fermion->odd, fermion->even,
                         &flops, &sent, &received, state->node);
    END_TIMING(state, flops, sent, received);

    return 0;
}    
