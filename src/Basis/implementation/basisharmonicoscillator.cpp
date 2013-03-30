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
    case 2:
        discretization2d();
        break;
    default:
        cerr << "BasisHydrogenLike::createInitalDiscretization():: "
             << "Dim = " << dim << " "
             << "not yet implemented." << endl;
        exit(EXIT_FAILURE);
    }
}
//------------------------------------------------------------------------------
void BasisHarmonicOscillator::discretization1d()
{
    C = field<cx_mat>(1);
    C(0) = zeros<cx_mat>(nGrid, nSpatialOrbitals);
    mat Ctmp(nGrid, nSpatialOrbitals);
    Wavefunction *wf;

    for(int j=0; j<nSpatialOrbitals; j++){
        wf = new HarmonicOscillator(cfg, states[2*j]);
        for(int i=0; i<nGrid; i++){
            Ctmp(i,j) = wf->evaluate(grid->x(i));
        }
        Ctmp.col(j) /= sqrt(dot(Ctmp.col(j),Ctmp.col(j)));

        delete wf;
    }

    C(0).set_real(Ctmp);

    // Forcing orthogonality by performing a SVD decomposition
    SVD(C(0));
#ifdef DEBUG
    cout << "BasisHarmonicOscillator::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            cout <<  "|C(" << i << ", " << j << ")| = " << sqrt(cdot(C(0).col(i),C(0).col(j))) << endl ;
        }
    }
#endif
}
//------------------------------------------------------------------------------
void BasisHarmonicOscillator::discretization2d()
{
    C = field<cx_mat>(nSpatialOrbitals);
    mat C2_tmp(nGrid, nGrid);

    for(int k=0; k<nSpatialOrbitals; k++){
        C(k) = zeros<cx_mat>(nGrid, nGrid);
        Wavefunction *wf = new HarmonicOscillator(cfg, states[2*k]);

        for(int i=0; i<nGrid; i++){
            for(int j=0; j<nGrid; j++){
                C2_tmp(i,j) = wf->evaluate(grid->x(i), grid->y(j));
            }
        }
        delete wf;

        C(k).set_real(C2_tmp);
        cx_mat Cnorm = conj(C(k)) % C(k);
        C(k) /= sqrt(accu(Cnorm));
    }

#ifdef DEBUG
    cout << "BasisHarmonicOscillator::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            cx_mat Cnorm = conj(C(i)) % C(j);
            cout <<  "psi_" << i << "* psi_" << j << " = " << accu(Cnorm) << endl ;
        }
    }
#endif
}
//------------------------------------------------------------------------------
