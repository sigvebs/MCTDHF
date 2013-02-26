#include "complextimerungekutta4.h"

//------------------------------------------------------------------------------
ComplexTimeRungeKutta4::ComplexTimeRungeKutta4(Config *cfg):
    ComplexTimeIntegrator(cfg)
{
//    i = cx_double(0,1);
    i = cx_double(1,0);
}
//------------------------------------------------------------------------------
void ComplexTimeRungeKutta4::stepForward()
{
    cout << "---------------------------------------------------------------\n";

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
    C = renormalize(C);
    A = A/sqrt(cdot(A,A));

    //--------------------------------------------------------------------------
    // Writing results to screen and file
    //--------------------------------------------------------------------------
    E = slater->getEnergy(A) ;
    if(std::isnan(E)){
        cerr << " E = Nan" << endl
             << "N = " << step << endl;
        exit(1);
    }

    // Saving C and A to disk
    EVec(step) = E;
    C.save("../DATA/C.mat", arma_ascii);
    A.save("../DATA/A.vec", arma_ascii);
    EVec.save("../DATA/EVec.mat", arma_ascii);

    cout.precision(16);
    cout << "step = " << step << endl
         << "E = " << E << endl
         << "|A|^2 = " << abs(cdot(A,A)) << endl
         << "|C(0)|^2 = " << abs(cdot(C.col(0),C.col(0))) << endl;

    for(uint i=1; i< C.n_cols; i++){
        cout << "|C("<<i<<")|^2 = " << abs(cdot(C.col(i),C.col(i))) << endl
             << "<C("<<i<<")| C(0)> = " << cdot(C.col(i),C.col(0)) << endl;
    }
    step++;
}
//------------------------------------------------------------------------------
cx_mat ComplexTimeRungeKutta4::renormalize(cx_mat C)
{
    // Re-normalization of C using SVD
    cx_mat X;
    vec s;
    cx_mat Y;
    svd_econ(X, s, Y, C);
    C = X*Y.t();

    return C;
}
//------------------------------------------------------------------------------
