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

    C = field<cx_mat>(1);
    C(0) = zeros<cx_mat>(nGrid, nSpatialOrbitals);

    mat Ctmp(nGrid, nSpatialOrbitals);
    Wavefunction* wf;

    for(int j=0; j<nSpatialOrbitals; j++){
        wf = new HydrogenLike(cfg, states[2*j]);
        for(int i=0; i<nGrid; i++){
            Ctmp(i,j) = wf->evaluate(grid->x(i));
        }
        Ctmp.col(j) /= sqrt(dot(Ctmp.col(j), Ctmp.col(j)));
        delete wf;
    }
    C(0).set_real(Ctmp);
    SVD(C(0));

#ifdef DEBUG
    cout << "BasisHydrogenLike::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            cout <<  "|C(" << i << ", " << j << ")| = " << sqrt(cdot(C(0).col(i), C(0).col(j))) << endl ;
        }
    }
#endif
}
//------------------------------------------------------------------------------
