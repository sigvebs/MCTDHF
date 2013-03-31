#include "randomunitarymatrix.h"

//------------------------------------------------------------------------------
RandomUnitaryMatrix::RandomUnitaryMatrix(Config *cfg):
    Basis(cfg)
{
}
//------------------------------------------------------------------------------
void RandomUnitaryMatrix::createInitalDiscretization()
{
    srand (time(NULL));
    C = randn<cx_mat>(nGrid, nSpatialOrbitals);
    SVD(C);
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
