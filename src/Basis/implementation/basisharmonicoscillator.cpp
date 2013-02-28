#include "basisharmonicoscillator.h"

//------------------------------------------------------------------------------
BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg):
    Basis(cfg)
{
    double latticeRange;
    try{
        nBasis = cfg->lookup("system.shells");
        latticeRange = cfg->lookup("spatialDiscretization.latticeRange");
        dx = cfg->lookup("spatialDiscretization.gridSpacing");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    nGrid = 2*latticeRange/dx+1;
    nSpatialOrbitals = states.size()/2;

    // Adding the number of gridpoints to the config file
    Setting &root = cfg->getRoot();
    Setting &tmp = root["spatialDiscretization"];
    tmp.add("nGrid", Setting::TypeInt) = nGrid;

    x = mat(nGrid, dim);
    for(int i=0; i<dim; i++)
        x.col(i) = linspace<vec>(-latticeRange,latticeRange,nGrid);
#ifdef DEBUG
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
        Ctmp.col(i) = wf->evaluate(x)*sqrt(dx);

        // Re-normalizing to remove numerical errors
        Ctmp.col(i) /= sqrt(dot(Ctmp.col(i),Ctmp.col(i)));
    }

    C.set_real(Ctmp);
    x.save(filnameAxis, arma_ascii);
    cout << "2" << endl;
#ifdef DEBUG
    cout << "BasisHarmonicOscillator::discretization1d()" << endl;
    for(int i=0; i<nSpatialOrbitals; i++){
        cout <<  "|C(" << i << ")| = " << sqrt(dot(Ctmp.col(i),Ctmp.col(i))) << endl ;
    }
#endif
}
//------------------------------------------------------------------------------
