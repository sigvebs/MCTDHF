#include "basisharmonicoscillator.h"

//------------------------------------------------------------------------------
BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg):
    Basis(cfg)
{
    double L;
    try{
        nBasis = cfg->lookup("system.shells");
        L = cfg->lookup("spatialDiscretization.latticeRange");
//        dx = cfg->lookup("spatialDiscretization.gridSpacing");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
//    nGrid = 2*latticeRange/dx+1;
    dx = 2.0*L/(double)(nGrid);
    nSpatialOrbitals = states.size()/2;

    // Adding the number of gridpoints to the config file
    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];

    tmp.add("gridSpacing", Setting::TypeFloat) = dx;
//    tmp.add("nGrid", Setting::TypeInt) = nGrid;

    x = mat(nGrid, dim);
    for(int i=0; i<dim; i++)
        x.col(i) = linspace<vec>(-L,L-dx,nGrid);

#ifdef DEBGUG
    cout << "BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg)" << endl
         << "nBasis \t\t= " << nBasis << endl
         << "latticeRange \t\t= " << latticeRange << endl
         << "dx \t\t= " << dx << endl
         << "nGrid \t\t= " << nGrid << endl;
#endif
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

        // Re-normalizing to remove numerical errors
        Ctmp.col(i) /= sqrt(dot(Ctmp.col(i),Ctmp.col(i)));
        delete wf;
    }

    C.set_real(Ctmp);
    x.save(filnameAxis, arma_ascii);
#ifdef DEBUG
    cout << "BasisHarmonicOscillator::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        cout <<  "|C(" << i << ")| = " << sqrt(dot(Ctmp.col(i),Ctmp.col(i))) << endl ;
    }
#endif
}
//------------------------------------------------------------------------------
