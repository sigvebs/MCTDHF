#include "complextimerungekutta4.h"

//------------------------------------------------------------------------------
ComplexTimeRungeKutta4::ComplexTimeRungeKutta4(Config *cfg):
    ComplexTimePropagation(cfg)
{
//    i = cx_double(0,1);
    i = cx_double(1,0);
}
//------------------------------------------------------------------------------
bool ComplexTimeRungeKutta4::stepForward()
{

#if 0
    // Testing -----------------
    if(step == 0)
        srand (time(NULL));
    cx_vec B = randu<cx_vec>( A.n_rows );
    cx_mat D = randn<cx_mat>(C.n_rows, C.n_cols);
    D.load("../DATA/C.mat");
    cx_mat X;
    vec s;
    cx_mat Y;
    svd_econ(X, s, Y, D);
    D = X*Y.t();

    cout << D << endl;
    B = B/sqrt(cdot(B, B));
    V->computeNewElements(D);
    h->computeNewElements(D);
    V->printInteractionElements();

    slater->computeRightHandSideComplexTime(B);
    cout << s << endl;
    cout <<"Done" << endl;
    exit(1);
    //-------------------
#else
    cx_vec k1, k2, k3, k4;
    cx_mat m1, m2, m3, m4;

    // Computing Runge-Kutta weights
    V->computeNewElements(C);
    h->computeNewElements(C);

    k1 = -i*dt*slater->computeRightHandSideComplexTime(A);
    m1 = -i*dt*orbital->computeRightHandSide(C, A);

    V->computeNewElements(C + 0.5*m1);
    h->computeNewElements(C + 0.5*m1);

    k2 = -i*dt*slater->computeRightHandSideComplexTime(A + 0.5*k1);
    m2 = -i*dt*orbital->computeRightHandSide(C + 0.5*m1, A + 0.5*k1);

    V->computeNewElements(C + 0.5*m2);
    h->computeNewElements(C + 0.5*m2);

    k3 = -i*dt*slater->computeRightHandSideComplexTime(A + 0.5*k2);
    m3 = -i*dt*orbital->computeRightHandSide(C + 0.5*m2, A + 0.5*k2);

    V->computeNewElements(C + m3);
    h->computeNewElements(C + m3);

    k4 = -i*dt*slater->computeRightHandSideComplexTime(A + k3);
    m4 = -i*dt*orbital->computeRightHandSide(C + m3, A + k3);

    // Computing new states
    A += 1.0/6.0*(k1 + 2*(k2 + k3) + k4);
    C += 1.0/6.0*(m1 + 2*(m2 + m3) + m4);
    t += dt;

    // Normalizing
    renormalize(C);
    A = A/sqrt(cdot(A,A));
#endif
    return true;
}
//------------------------------------------------------------------------------
