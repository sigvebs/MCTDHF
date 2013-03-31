#include "basishydrogenlike.h"
//------------------------------------------------------------------------------
BasisHydrogenLike::BasisHydrogenLike(Config *cfg):
    Basis(cfg)
{
}
//------------------------------------------------------------------------------
void BasisHydrogenLike::createInitalDiscretization()
{
    C = zeros<cx_mat>(nGrid, nSpatialOrbitals);

    mat Ctmp(nGrid, nSpatialOrbitals);
    Wavefunction* wf;

    for(int j=0; j<nSpatialOrbitals; j++){
        wf = new HydrogenLike(cfg, states[2*j]);
        for(int i=0; i<nGrid; i++){
            Ctmp(i,j) = wf->evaluate(grid->at(i));
        }
        delete wf;
    }
    C.set_real(Ctmp);
    SVD(C);

#ifdef DEBUG
    cout << "BasisHydrogenLike::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            cout <<  "|C(" << i << ", " << j << ")| = " << sqrt(cdot(C.col(i), C.col(j))) << endl ;
        }
    }
#endif
}
//------------------------------------------------------------------------------
