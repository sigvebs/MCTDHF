#include "complextimecranknicholson.h"
//------------------------------------------------------------------------------
ComplexTimeCrankNicholson::ComplexTimeCrankNicholson(Config *cfg):
    ComplexTimePropagation(cfg)
{
}
//------------------------------------------------------------------------------
bool ComplexTimeCrankNicholson::stepForward()
{
    /*
    n = H.n_cols;
    I = eye<cx_mat>(n,n);

    H1 = I + 0.5*dt*H;
    H2 = I - 0.5*dt*H;

    H1 = inv(H1);

    A = H1*H2*A;
    t += dt;
    */
    return 1;
}
//------------------------------------------------------------------------------
