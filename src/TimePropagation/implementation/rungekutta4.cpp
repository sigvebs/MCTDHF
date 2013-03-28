#include "rungekutta4.h"

//------------------------------------------------------------------------------
RungeKutta4::RungeKutta4(Config *cfg):
    TimePropagation(cfg)
{
    i = cx_double(0,1);
}
//------------------------------------------------------------------------------
bool RungeKutta4::stepForward()
{
    // Computing Runge-Kutta weights
    V->computeNewElements(C);
    h->computeNewElements(C, t);

    k1 = -i*dt*slater->computeRightHandSide(A);
    m1 = -i*dt*orbital->computeRightHandSide(C, A);

    V->computeNewElements(C + 0.5*m1);
    h->computeNewElements(C + 0.5*m1, t);

    k2 = -i*dt*slater->computeRightHandSide(A + 0.5*k1);
    m2 = -i*dt*orbital->computeRightHandSide(C + 0.5*m1, A + 0.5*k1);

    V->computeNewElements(C + 0.5*m2);
    h->computeNewElements(C + 0.5*m2, t);

    k3 = -i*dt*slater->computeRightHandSide(A + 0.5*k2);
    m3 = -i*dt*orbital->computeRightHandSide(C + 0.5*m2, A + 0.5*k2);

    V->computeNewElements(C + m3);
    h->computeNewElements(C + m3, t);

    k4 = -i*dt*slater->computeRightHandSide(A + k3);
    m4 = -i*dt*orbital->computeRightHandSide(C + m3, A + k3);

    // Computing new states
    A += 1.0/6.0*(k1 + 2*(k2 + k3) + k4);
    C += 1.0/6.0*(m1 + 2*(m2 + m3) + m4);

    // Normalizing
    renormalize(C);
    A = A/sqrt(cdot(A,A));

    return 1;
}
//------------------------------------------------------------------------------
