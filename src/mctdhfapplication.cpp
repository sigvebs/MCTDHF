#include "mctdhfapplication.h"

//------------------------------------------------------------------------------
MctdhfApplication::MctdhfApplication(string configFilename)
{
    cout << "Loading " << configFilename << endl;

    try {
        cfg.readFile(configFilename.c_str());
    } catch (const FileIOException &fioex) {
        cerr << "I/O error while reading config file." << endl;
        exit(EXIT_FAILURE);
    }

    try {
        loadDataset = cfg.lookup("systemSettings.loadDataset");
        cfg.lookupValue("systemSettings.loadDatasetPath", loadDatasetPath);
    }catch (const SettingNotFoundException &fioex) {
        loadDataset = false;
        loadDatasetPath = "";
    }

#ifdef DEBUG
    cout << "Reading configuration from file" << endl;
#endif
}
//------------------------------------------------------------------------------
void MctdhfApplication::run()
{

    //--------------------------------------------------------------------------
    // Creating the orbitals and spatial discretization
    //--------------------------------------------------------------------------
    Basis* orb = setBasis();
    orb->createBasis();
    const vector<vec> orbitals = orb->getBasis();
    vec x;
    cx_mat C;
    if(loadDataset){
        x.load(loadDatasetPath + "/x.mat");
        C.load(loadDatasetPath + "/C.mat");

        checkXC(x, C, orbitals);
    }else{
        orb->createInitalDiscretization();
        C = orb->getInitalOrbitals();
        x = orb->getX();
    }
    delete orb;

    //--------------------------------------------------------------------------
    // Creating all possible Slater determinants from the set of orbitals
    //--------------------------------------------------------------------------
    SlaterDeterminants slater(&cfg, orbitals);
    slater.createSlaterDeterminants();
    const vector<bitset<BITS> > &slaterDeterminants = slater.getSlaterDeterminants();

    //--------------------------------------------------------------------------
    // Interaction operator
    //--------------------------------------------------------------------------
    cout << "Setting up the interaction operator" << endl;
    Interaction V(&cfg, setMeanFieldIntegrator());
    setInteractionPotentials(V, x);

    //--------------------------------------------------------------------------
    // Setting the single particle operator
    //--------------------------------------------------------------------------
    cout << "Setting up the single particle operator" << endl;
    SingleParticleOperator h(&cfg, setDifferentialOpertor(x));
    setOneBodyPotentials(h, x);

    //--------------------------------------------------------------------------
    // Setting up the Slater equation
    //--------------------------------------------------------------------------
    cout << "Setting up the Slater determiant equation" << endl;
    SlaterEquation slaterEquation(&cfg, slaterDeterminants, &V, &h);

    //--------------------------------------------------------------------------
    // Loading/Creating an initial coefficient vector for the Slater determinants.
    //--------------------------------------------------------------------------
    cx_vec A;
    if(loadDataset){
        A.load(loadDatasetPath + "/A.mat");
    }else{
        A = randu<cx_vec>(slaterDeterminants.size());
        A = A/sqrt(cdot(A, A));
    }

    //--------------------------------------------------------------------------
    // Setting up the orbital equation
    //--------------------------------------------------------------------------
    cout << "Setting up the Orbital equation" << endl;
    OrbitalEquation orbEq(&cfg, slaterDeterminants, &V, &h);

    //--------------------------------------------------------------------------
    // Setting up the complex time integrator
    //--------------------------------------------------------------------------
    ComplexTimePropagation* complexTimePropagation = setComplexTimeIntegrator();
    complexTimePropagation->setDependencies(&slaterEquation, &orbEq, &V, &h);
    complexTimePropagation->setInititalState(A, C);

    //--------------------------------------------------------------------------
    // Performing complex time integrator
    //--------------------------------------------------------------------------
    cout << "Starting imaginary time propagation" << endl;
    complexTimePropagation->doComplexTimePropagation();
    A = complexTimePropagation->getCurrentA();
    C = complexTimePropagation->getCurrentC();
    delete complexTimePropagation;

    //--------------------------------------------------------------------------
    // Time integration
    //--------------------------------------------------------------------------
    bool doTimeIntegration = cfg.lookup("systemSettings.doTimeIntegration");
    if(doTimeIntegration){
        TimePropagation *timePropagator = setTimeIntegrator();
        setTimeDepOneBodyPotentials(h, x);
        timePropagator->setDependencies(&slaterEquation, &orbEq, &V, &h);
        timePropagator->setInititalState(A, C);

        cout << "Starting time propagation" << endl;
        timePropagator->doComplexTimePropagation();

        delete timePropagator;
    }
    //--------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
void MctdhfApplication::setInteractionPotentials(Interaction &V, const vec &x)
{
    InteractionPotential* I;
    int interactionType = cfg.lookup("interactionPotential.interactionType");

    switch (interactionType) {
    case IP_HARMONIC_OSCILLATOR:
        I = new HarmonicOscillatorInteraction(&cfg, x);
        V.addPotential(I);
        break;
    case IP_SHEILDED_COULOMB:
        I = new ScreenedCoulombInteraction(&cfg, x);
        V.addPotential(I);
        break;
    default:
        cerr << "Interaction not implemented:: " << interactionType << endl;
        exit(EXIT_FAILURE);
    }

    V.updatePositionBasisElements();
}
//------------------------------------------------------------------------------
void MctdhfApplication::setOneBodyPotentials(SingleParticleOperator &h, const vec &x)
{
    Potential* I;
    const Setting& root = cfg.getRoot();
    const Setting &oneBodyPotentials = root["oneBodyPotential"]["potential"];
    int nPotentials = oneBodyPotentials.getLength();

    for(int i=0; i<nPotentials; i++){
        int potential = oneBodyPotentials[i];
        switch (potential) {
        case HARMONIC_OSCILLATOR_ONE_BODY:
            I = new HarmonicOscillatorOneBody(&cfg, x);
            h.addPotential(I);
            break;
        case COULOMB_INTERACTION_NUCLEUS:
            I = new CoulombInteractionNucleus(&cfg, x);
            h.addPotential(I);
            break;
        case ANHARMONIC_DOUBLE_WELL:
            I = new AnharmonicDoubleWell(&cfg, x);
            h.addPotential(I);
            break;
        default:
            cerr << "Potential not implemented:: " << potential << endl;
            exit(EXIT_FAILURE);
        }
    }
}
//------------------------------------------------------------------------------
void MctdhfApplication::setTimeDepOneBodyPotentials(SingleParticleOperator &h, const vec &x)
{
    Potential* I;

    const Setting& root = cfg.getRoot();
    const Setting &oneBodyPotentials = root["oneBodyPotential"]["timeDepPotential"];
    int nPotentials = oneBodyPotentials.getLength();

    for(int i=0; i<nPotentials; i++){
        int potential = oneBodyPotentials[i];

        switch (potential) {
        case SIMPLE_LASER:
            I = new simpleLaser(&cfg, x);
            h.addPotential(I);
            break;
        default:
            cerr << "Time dependent potential not implemented:: " << potential << endl;
            exit(EXIT_FAILURE);
        }
    }
}
//------------------------------------------------------------------------------
TimePropagation *MctdhfApplication::setTimeIntegrator()
{
    // Setting the integrator
    TimePropagation *I;
    int integrator = cfg.lookup("timeIntegration.integrator");

    switch (integrator) {
    case RUNGE_KUTTA4:
        I = new RungeKutta4(&cfg);
        break;
    case RUNGE_KUTTA_FEHLBERG :
        I = new RungeKutta4(&cfg);
        break;
    default:
        cerr << "Time Integrator not implemented:: " << integrator << endl;
        exit(EXIT_FAILURE);
    }

    return I;
}
//------------------------------------------------------------------------------
ComplexTimePropagation* MctdhfApplication::setComplexTimeIntegrator()
{
    // Setting the integrator
    ComplexTimePropagation *I;
    int integrator = cfg.lookup("ComplexTimeIntegration.integrator");

    switch (integrator) {
    case CT_CRANK_NICOLSON:
        I = new ComplexTimeCrankNicholson(&cfg);
        break;
    case CT_RUNGE_KUTTA4:
        I = new ComplexTimeRungeKutta4(&cfg);
        break;
    case CT_RUNGE_KUTTA_FEHLBERG :
        I = new ComplexTimeRungeKuttaFehlberg(&cfg);
        break;
    default:
        cerr << "Complex Time Integrator not implemented:: " << integrator << endl;
        exit(EXIT_FAILURE);
    }

    return I;
}
//------------------------------------------------------------------------------
DifferentialOperator* MctdhfApplication::setDifferentialOpertor(const vec &x)
{
    // Setting the integrator
    DifferentialOperator *I;
    int differentialOperator = cfg.lookup("spatialDiscretization.differentialOperator");

    switch (differentialOperator) {
    case DO_FINITE_DIFFERENCE_1d:
        I = new FiniteDifference1d(&cfg, x);
        break;
    case DO_FINITE_DIFFERENCE_FIVE_POINT_1D:
        I = new FiniteDifferenceFivePoint1d(&cfg, x);
        break;
    case DO_SPECTRAL_1D:
        I = new Spectral1d(&cfg, x);
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
    case OBT_RAND_UNITARY_MATRIX:
        I = new RandomUnitaryMatrix(&cfg);
        break;
    default:
        cerr << "Basis not implemented:: " << basisType << endl;
        exit(EXIT_FAILURE);
    }
    return I;
}
//------------------------------------------------------------------------------
void MctdhfApplication::checkXC(const vec &x, const cx_mat &C, const vector<vec> &orbitals)
{
    // Checking that x has the right dimensions
    int nGrid = cfg.lookup("spatialDiscretization.nGrid");
    if(x.n_elem != nGrid){
        cerr << "The loaded x-matrix does not have the right dimensions" << endl;
        cerr << "elements(x): " << x.n_elem << endl;
        cerr << "config nGrid: = " << nGrid << endl;
        exit(EXIT_FAILURE);
    }

    // Checking that C has the right dimensions
    if(C.n_cols != orbitals.size()/2){
        cerr << "The loaded C-matrix does not have the right dimensions" << endl;
        cerr << "Orbitals(C): " << 2*C.n_cols << endl;
        cerr << "orbitlas config: = " << orbitals.size()<< endl;
        exit(EXIT_FAILURE);
    }
}
//-----------------------------------------------------------------------------
