#include "randomunitarymatrix.h"

//------------------------------------------------------------------------------
RandomUnitaryMatrix::RandomUnitaryMatrix(Config *cfg):
    Basis(cfg)
{
}
//------------------------------------------------------------------------------
void RandomUnitaryMatrix::createInitalDiscretization()
{
    switch(dim)
    {
    case 1:
        discretization1d();
        break;
    default:
        cerr << "RandomUnitaryMatrix::createInitalDiscretization():: "
             << "Dim = " << dim << " "
             << "not yet implemented." << endl;
        exit(EXIT_FAILURE);
    }
}
//------------------------------------------------------------------------------
void RandomUnitaryMatrix::discretization1d()
{
    srand (time(NULL));
    C = randn<cx_mat>(nGrid, nSpatialOrbitals);
    cx_mat X;
    vec s;
    cx_mat Y;
    svd_econ(X, s, Y, C);
    C = X*Y.t();
#ifdef DEBUG
    cout << "RandomUnitaryMatrix::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            cout <<  "|C(" << i << ", " << j << ")| = " << sqrt(cdot(C.col(i),C.col(j))) << endl ;
        }
    }
#endif
}
//------------------------------------------------------------------------------
