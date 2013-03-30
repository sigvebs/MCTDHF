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

    C = field<cx_mat>(1);
    C(0) = randn<cx_mat>(nGrid, nSpatialOrbitals);
    SVD(C(0));
#ifdef DEBUG
    cout << "RandomUnitaryMatrix::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            cout <<  "|C(" << i << ", " << j << ")| = " << sqrt(cdot(C(0).col(i),C(0).col(j))) << endl ;
        }
    }
#endif
}
//------------------------------------------------------------------------------
