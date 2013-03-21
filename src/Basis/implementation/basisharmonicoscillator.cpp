#include "basisharmonicoscillator.h"

//------------------------------------------------------------------------------
BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg):
    Basis(cfg)
{
}
//------------------------------------------------------------------------------
void BasisHarmonicOscillator::createInitalDiscretization()
{
    switch(dim)
    {
    case 1:
        discretization1d();
        break;
    default:
        cerr << "BasisHarmonicOscillator::createInitalDiscretization():: "
             << "Dim = " << dim << " "
             << "not yet implemented." << endl;
        exit(EXIT_FAILURE);
    }
}
//------------------------------------------------------------------------------
void BasisHarmonicOscillator::discretization1d()
{
    cout << "BasisHarmonicOscillator" << endl;
    C = zeros<cx_mat>(nGrid, nSpatialOrbitals);
    mat Ctmp(nGrid, nSpatialOrbitals);
    Wavefunction* wf;

    for(int i=0; i<nSpatialOrbitals; i++){
        wf = new HarmonicOscillator1d(cfg, states[2*i]);
        Ctmp.col(i) = wf->evaluate(x);
        Ctmp.col(i) /= sqrt(dot(Ctmp.col(i),Ctmp.col(i)));
        delete wf;
    }

    C.set_real(Ctmp);

    // Forcing orthogonality by performing a SVD decomposition
    cx_mat X;
    vec s;
    cx_mat Y;
    svd_econ(X, s, Y, C);
    C = X*Y.t();
#ifdef DEBUG
    cout << "BasisHarmonicOscillator::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            cout <<  "|C(" << i << ", " << j << ")| = " << sqrt(cdot(C.col(i),C.col(j))) << endl ;
        }
    }
#endif
}
//------------------------------------------------------------------------------
