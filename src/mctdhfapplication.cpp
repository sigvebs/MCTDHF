#include "mctdhfapplication.h"

//------------------------------------------------------------------------------
MctdhfApplication::MctdhfApplication(int* argc, char ***argv, string configFilename)
{
    cout << "Loading " << configFilename << endl;

    try {
        cfg.readFile(configFilename.c_str());
    } catch (const FileIOException &fioex) {
        cerr << "I/O error while reading config file." << endl;
        exit(EXIT_FAILURE);
    }

#ifdef DEBUG
    cout << "Reading configuration from file" << endl;
#endif
}
//------------------------------------------------------------------------------
void MctdhfApplication::run()
{
    // Creating the orbitals
    Basis* orb = new BasisHarmonicOscillator(&cfg);
    orb->createBasis();
    orb->createInitalDiscretization();
    const vector<vec> &orbitals = orb->getBasis();
    const cx_mat &C = orb->getInitalOrbitals();

    // Creating all possible Slater determinants from the set of orbitals
    SlaterDeterminants slater(&cfg, orbitals);
    slater.createSlaterDeterminants();
    const vector<bitset<BITS> > &slaterDeterminants = slater.getSlaterDeterminants();

    // Interaction operator
    cout << "Setting up the interaction operator" << endl;
    Interaction V(&cfg, orbitals);
    V.computeNewElements(C);

    // Setting the single particle operator
    cout << "Setting up the single particle operator" << endl;
    SingleParticleOperator h(&cfg, orbitals);
    h.computeNewElements(C);

    // Setting up the Slater equation
    cout << "Setting up the Slater determiant equation" << endl;
    SlaterEquation slaterEquation(&cfg, orbitals, slaterDeterminants, &V, &h);
    slaterEquation.setInitalState();
    cout << "E = " << slaterEquation.getEnergy() << endl;

    // Setting up the orbital equation
    cout << "Setting up the Orbital equation" << endl;
    OrbitalEquation orbEq(&cfg, orbitals, slaterDeterminants, &V, &h);
    orbEq.setInititalState(C);


    // Setting up the complex time integrator
    ComplexTimeIntegrator* complexTimeIntegrator = setComplexTimeIntegrator();
    complexTimeIntegrator->setDependencies(&slaterEquation, &orbEq, &V, &h);
    ComplexTimePropagation complexTimePropagation(&cfg, complexTimeIntegrator);

    cout << "Complex time propagation" << endl;
    complexTimePropagation.doComplexTimePropagation();

    // Cleaning memory
    delete orb;

    cout << "Done" << endl;
}
//------------------------------------------------------------------------------
ComplexTimeIntegrator* MctdhfApplication::setComplexTimeIntegrator()
{
    // Setting the integrator
    ComplexTimeIntegrator *I;
    int complexTimeIntegrator = cfg.lookup("ComplexTimeIntegration.integrator");

    switch (complexTimeIntegrator) {
    case CT_CRANK_NICOLSON:
        I = new ComplexTimeCrankNicholson(&cfg);
        break;
    case CT_RUNGE_KUTTA4:
        I = new ComplexTimeRungeKutta4(&cfg);
        break;
    }

    return I;
}
//-----------------------------------------------------------------------------
//WaveFunction* MctdhfApplication::setWavefunction()
//{
//    WaveFunction *wf;
//    int dim = cfg.lookup("system.dim");

//    switch(dim){
//    case 1:
//        wf = new HarmonicOscillator1d(&cfg);
//        break;
//    case 2:
//        wf = new HarmonicOscillator2d(&cfg);
//        break;
//    }

//    return wf;
//}
////------------------------------------------------------------------------------
//SpatialIntegrator* MctdhfApplication::setSpatialIntegrator(WaveFunction* wf)
//{
//    // Setting the integrator
//    SpatialIntegrator *I;
//    int spatialIntegrator;
//    cfg.lookupValue("spatialIntegration.integrator", spatialIntegrator);

//    switch (spatialIntegrator) {
//    case MONTE_CARLO:
//        I = new MonteCarloIntegrator(&cfg);
//        break;
//    case GAUSS_LAGUERRE:
//        I = new GaussLaguerreIntegrator(&cfg, wf);
//        break;
//    case GAUSS_HERMITE:
//        I = new GaussHermiteIntegrator(&cfg);
//        break;
//    case INTERACTION_INTEGRATOR:
//        I = new InteractonIntegrator(&cfg);
//        break;
//    case MONTE_CARLO_IS:
//        I = new MonteCarloImportanceSampled(&cfg);
//        break;
//    }

//    return I;
//}
////------------------------------------------------------------------------------

