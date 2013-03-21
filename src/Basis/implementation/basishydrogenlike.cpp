#include "basishydrogenlike.h"
//------------------------------------------------------------------------------
BasisHydrogenLike::BasisHydrogenLike(Config *cfg):
    Basis(cfg)
{
}
//------------------------------------------------------------------------------
void BasisHydrogenLike::createInitalDiscretization()
{
    switch(dim)
    {
    case 1:
        discretization1d();
        break;
    default:
        cerr << "BasisHydrogenLike::createInitalDiscretization():: "
             << "Dim = " << dim << " "
             << "not yet implemented." << endl;
        exit(EXIT_FAILURE);
    }
}
//------------------------------------------------------------------------------
void BasisHydrogenLike::discretization1d()
{ cout << "discreet" << endl;
    C = zeros<cx_mat>(nGrid, nSpatialOrbitals);
    mat Ctmp(nGrid, nSpatialOrbitals);
    Wavefunction* wf;

    for(int i=0; i<nSpatialOrbitals; i++){
        wf = new HydrogenLike(cfg, states[2*i]);
        Ctmp.col(i) = wf->evaluate(x);
    }

    C.set_real(Ctmp);

    // Forcing orthogonality by performing a SVD decomposition
    cx_mat X;
    vec s;
    cx_mat Y;
    svd_econ(X, s, Y, C);
    C = X*Y.t();

#ifdef DEBUG
    cout << "BasisHydrogenLike::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            cout <<  "|C(" << i << ", " << j << ")| = " << sqrt(cdot(C.col(i),C.col(j))) << endl ;
        }
    }
#endif
}
//------------------------------------------------------------------------------
