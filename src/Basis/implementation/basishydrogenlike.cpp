#include "basishydrogenlike.h"
//------------------------------------------------------------------------------
BasisHydrogenLike::BasisHydrogenLike(Config *cfg):
    Basis(cfg)
{
    double L;
    try{
        nBasis = cfg->lookup("system.shells");
        L = cfg->lookup("spatialDiscretization.latticeRange");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    dx = 2.0*L/(double)(nGrid);
    nSpatialOrbitals = states.size()/2;

    // Adding the number of gridpoints to the config file
    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];

    tmp.add("gridSpacing", Setting::TypeFloat) = dx;

    x = mat(nGrid, dim);
    for(int i=0; i<dim; i++)
        x.col(i) = linspace<vec>(-L,L-dx,nGrid);

#ifdef DEBUG
    cout << "BasisHydrogenLike::BasisHydrogenLike(Config *cfg)" << endl
         << "nBasis \t\t= " << nBasis << endl
         << "latticeRange \t\t= " << latticeRange << endl
         << "dx \t\t= " << dx << endl
         << "nGrid \t\t= " << nGrid << endl;
#endif
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

        // Re-normalizing to remove numerical errors
        Ctmp.col(i) = Ctmp.col(i)/sqrt(dot(Ctmp.col(i),Ctmp.col(i)));
    }

    C.set_real(Ctmp);
    x.save(filnameAxis, arma_ascii);

#ifdef DEBUG
//#if 1
    cout << "BasisHarmonicOscillator::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        cout <<  "|C(" << i << ")|^2 = " << cdot(C.col(i),C.col(i)) << endl ;
    }
#endif
}
//------------------------------------------------------------------------------
