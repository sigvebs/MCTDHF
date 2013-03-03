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
    Basis* orb = setBasis();
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
    Interaction V(&cfg, setMeanFieldIntegrator());
    setInteractionPotentials(V);

    // Setting the single particle operator
    cout << "Setting up the single particle operator" << endl;
    SingleParticleOperator h(&cfg, setDifferentialOpertor());

    // Setting up the Slater equation
    cout << "Setting up the Slater determiant equation" << endl;
    SlaterEquation slaterEquation(&cfg, slaterDeterminants, &V, &h);

    // Creating an initial coefficient vector for the Slater determinants.
    cx_vec A = randu<cx_vec>(slaterDeterminants.size());
    A = A/sqrt(cdot(A, A));

    // Setting up the orbital equation
    cout << "Setting up the Orbital equation" << endl;
    OrbitalEquation orbEq(&cfg, slaterDeterminants, &V, &h);

    // Setting up the complex time integrator
    ComplexTimeIntegrator* complexTimeIntegrator = setComplexTimeIntegrator();
    complexTimeIntegrator->setDependencies(&slaterEquation, &orbEq, &V, &h);
    complexTimeIntegrator->setInititalState(A, C);
    ComplexTimePropagation complexTimePropagation(&cfg, complexTimeIntegrator);

    cout << "Complex time propagation" << endl;
    complexTimePropagation.doComplexTimePropagation();

    // Cleaning memory
    delete orb;
    delete complexTimeIntegrator;

    cout << "Done" << endl;
}
//------------------------------------------------------------------------------
void MctdhfApplication::setInteractionPotentials(Interaction &V)
{
    InteractionPotential* I;
    int interactionType = cfg.lookup("interactionPotential.interactionType");

    switch (interactionType) {
    case IP_HARMONIC_OSCILLATOR:
        I = new HarmonicOscillatorInteraction(&cfg);
        V.addPotential(I);
        break;
    case IP_SHEILDED_COULOMB:
        I = new ScreenedCoulombInteraction(&cfg);
        V.addPotential(I);
        break;
    default:
        cerr << "Interaction not implemented:: " << interactionType << endl;
        exit(EXIT_FAILURE);
    }

    V.updatePositionBasisElements();
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
    default:
        cerr << "Complex Time Integrator not implemented:: " << complexTimeIntegrator << endl;
        exit(EXIT_FAILURE);
    }

    return I;
}
//------------------------------------------------------------------------------
DifferentialOperator* MctdhfApplication::setDifferentialOpertor()
{
    // Setting the integrator
    DifferentialOperator *I;
    int differentialOperator = cfg.lookup("spatialDiscretization.differentialOperator");

    switch (differentialOperator) {
    case DO_FINITE_DIFFERENCE_1d:
        I = new FiniteDifference1d(&cfg);
        break;
    case DO_FINITE_DIFFERENCE_FIVE_POINT_1D:
        I = new FiniteDifferenceFivePoint1d(&cfg);
        break;
    case DO_SPECTRAL_1D:
        I = new Spectral1d(&cfg);
        break;
    default:
        cerr << "Differential Operator not implemented:: " << differentialOperator << endl;
        exit(EXIT_FAILURE);
    }
    return I;
}
//------------------------------------------------------------------------------
MeanFieldIntegrator *MctdhfApplication::setMeanFieldIntegrator()
{
    MeanFieldIntegrator* I;
    int meanFieldIntegrator = cfg.lookup("meanFieldIntegrator.integratorType");

    switch (meanFieldIntegrator) {
    case MF_TRAPEZODIAL:
        I = new MfTrapezoidal(&cfg);
        break;
    case MF_LOW_RANK_APPROXIMATION:
        I = new MfLowRankApproximation(&cfg);
        break;
    default:
        cerr << "Mean Field integrator not implemented:: " << meanFieldIntegrator << endl;
        exit(EXIT_FAILURE);
    }
    return I;
}
//------------------------------------------------------------------------------
Basis *MctdhfApplication::setBasis()
{
    // Setting the basis
    Basis *I;
    int basisType = cfg.lookup("wavefunction.basisType");

    switch (basisType) {
    case OBT_HARMONIC_OSCILLATOR:
        I = new BasisHarmonicOscillator(&cfg);
        break;
    case OBT_HYDROGEN_LIKE:
        I = new BasisHydrogenLike(&cfg);
        break;
    default:
        cerr << "Basis not implemented:: " << basisType << endl;
        exit(EXIT_FAILURE);
    }
    return I;
}
//-----------------------------------------------------------------------------
