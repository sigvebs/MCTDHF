#include "complextimerungekutta4.h"

//------------------------------------------------------------------------------
ComplexTimeRungeKutta4::ComplexTimeRungeKutta4(Config *cfg):
    ComplexTimeIntegrator(cfg)
{
    Eprev = E = 0;

    // Tmp
    int N = cfg->lookup("ComplexTimeIntegration.N");
    EVec = zeros(N);
    step = 0;
    i = cx_double(1,0);
}
//------------------------------------------------------------------------------
void ComplexTimeRungeKutta4::stepForward()
{
    cout << "-------------------------------------\n";

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

    double E = slater->getEnergy() ;
    Eprev = E;
    if(std::isnan(E)){
        cerr << " E = Nan" << endl
             << "N = " << step << endl;
        exit(1);
    }


    // Saving C and A to disk
    mat Ctmp = real(C);
    EVec(step) = E;
    Ctmp.save("../DATA/C.mat", arma_ascii);
    A.save("../DATA/A.vec", arma_ascii);
    EVec.save("../DATA/EVec.mat", arma_ascii);

//    C = renomralize(C);
    int nSpatialOrbitals = C.n_cols;
    cx_mat O = zeros<cx_mat>(nSpatialOrbitals,nSpatialOrbitals);
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int l=i; l<nSpatialOrbitals; l++){
            O(i,l) = cdot(C.col(i), C.col(l));
            O(l,i) = conj(O(i,l));
        }
    }
    cout << O << endl;

    cout << "E = " << E << "\t\t step = " << step << endl
         << "norm(A) = " << cdot(A,A) << endl;



    step++;
}
//------------------------------------------------------------------------------
cx_mat ComplexTimeRungeKutta4::renomralize(cx_mat C)
{
    for(int i=0; i<C.n_cols; i++){
        C.col(i) /= cdot(C.col(i), C.col(i));
    }

    return C;
}
//------------------------------------------------------------------------------
