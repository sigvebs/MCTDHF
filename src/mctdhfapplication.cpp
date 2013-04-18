#include "mctdhfapplication.h"

//------------------------------------------------------------------------------
MctdhfApplication::MctdhfApplication(string configFilename)
{
    //--------------------------------------------------------------------------
    // Setting up MPI
    //--------------------------------------------------------------------------
    myRank = 0;
#ifdef USE_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif
    isMaster = (bool)(myRank == 0);

    cout << "Loading " << configFilename << endl;

    try {
        cfg.readFile(configFilename.c_str());
    } catch (const FileIOException &fioex) {
        cerr << "I/O error while reading config file." << endl;
        exit(EXIT_FAILURE);
    }

    try {
        loadDataset = cfg.lookup("systemSettings.loadDataset");
        doTimeIntegration = cfg.lookup("systemSettings.doTimeIntegration");
        doComplexTimeIntegration = cfg.lookup("systemSettings.doComplexTimeIntegration");
        cfg.lookupValue("systemSettings.loadDatasetPath", loadDatasetPath);
    }catch (const SettingNotFoundException &fioex) {
        loadDataset = false;
    }

#ifdef DEBUG
    cout << "Reading configuration from file" << endl;
#endif
}
//------------------------------------------------------------------------------
void MctdhfApplication::run()
{
    //--------------------------------------------------------------------------
    // Setting up the spatial discretization
    //--------------------------------------------------------------------------
    Grid grid(&cfg);
    if(loadDataset){
        grid.loadGrid();
    }else{
        grid.createInitalDiscretization();
    }
    if(isMaster)
        grid.saveGrid();

    //--------------------------------------------------------------------------
    // Setting up the orbitals
    //--------------------------------------------------------------------------
    Basis *orbitalBasis;
    cx_mat C;

    orbitalBasis = setBasis();
    orbitalBasis->setGrid(&grid);

    if(loadDataset){
        orbitalBasis->loadOrbitals();
    }else{
        orbitalBasis->createBasis();
        orbitalBasis->createInitalDiscretization();
    }
    const vector<vec> orbitals = orbitalBasis->getBasis();

    C = orbitalBasis->getOrbitals();
    delete orbitalBasis;

    //--------------------------------------------------------------------------
    // Creating all possible Slater determinants from the set of orbitals
    //--------------------------------------------------------------------------
    SlaterDeterminants slater(&cfg, orbitals);

    if(loadDataset){
        slater.load();
    }else{
        slater.createSlaterDeterminants();
        slater.createInitialState();
    }
    const vector<bitset<BITS> > &slaterDeterminants = slater.getSlaterDeterminants();
    cx_vec A = slater.getCoefficients();

    //--------------------------------------------------------------------------
    // Interaction operator
    //--------------------------------------------------------------------------
    if(isMaster)
        cout << "Setting up the interaction operator" << endl;
    Interaction V(&cfg, setMeanFieldIntegrator());
    setInteractionPotentials(V, grid);

    //--------------------------------------------------------------------------
    // Setting the single particle operator
    //--------------------------------------------------------------------------
    if(isMaster)
        cout << "Setting up the single particle operator" << endl;
    SingleParticleOperator h(&cfg, setDifferentialOpertor(grid));
    setOneBodyPotentials(h, grid);

    //--------------------------------------------------------------------------
    // Setting up the Slater equation
    //--------------------------------------------------------------------------
    if(isMaster)
        cout << "Setting up the Slater determiant equation" << endl;
    SlaterEquation slaterEquation(&cfg, slaterDeterminants, &V, &h);

    //--------------------------------------------------------------------------
    // Setting up the orbital equation
    //--------------------------------------------------------------------------
    if(isMaster)
        cout << "Setting up the Orbital equation" << endl;
    OrbitalEquation orbEq(&cfg, slaterDeterminants, &V, &h);

    //--------------------------------------------------------------------------
    // Setting up and performing complex time integration
    //--------------------------------------------------------------------------
    if(doComplexTimeIntegration){
        ComplexTimePropagation* complexTimePropagation = setComplexTimeIntegrator();
        complexTimePropagation->setDependencies(&slaterEquation, &orbEq, &V, &h);
        complexTimePropagation->setInititalState(A, C);

        cout << "Starting imaginary time propagation" << endl;
        complexTimePropagation->doComplexTimePropagation();
        A = complexTimePropagation->getCurrentA();
        C = complexTimePropagation->getCurrentC();
        delete complexTimePropagation;
    }

    //--------------------------------------------------------------------------
    // Time integration
    //--------------------------------------------------------------------------
    if(doTimeIntegration){
        TimePropagation *timePropagator = setTimeIntegrator();
        setTimeDepOneBodyPotentials(h, grid);
        timePropagator->setDependencies(&slaterEquation, &orbEq, &V, &h);
        timePropagator->setInititalState(A, C);

        cout << "Starting time propagation" << endl;
        timePropagator->doTimePropagation();
        delete timePropagator;
    }
    //--------------------------------------------------------------------------
#ifdef USE_MPI
    MPI_Finalize();
#endif
}
//------------------------------------------------------------------------------
void MctdhfApplication::setInteractionPotentials(Interaction &V, const Grid &grid)
{
    InteractionPotential* I;
    const Setting& root = cfg.getRoot();
    const Setting &interactionPotentials = root["interactionPotential"]["interactionType"];
    int nPotentials = interactionPotentials.getLength();
//    int interactionType = cfg.lookup("interactionPotential.interactionType");

    for(int i=0; i<nPotentials; i++){
        int interactionType = interactionPotentials[i];
        switch (interactionType) {
        case IP_HARMONIC_OSCILLATOR:
            I = new HarmonicOscillatorInteraction(&cfg, grid);
            V.addPotential(I);
            break;
        case IP_SHIELDED_COULOMB:
            I = new ScreenedCoulombInteraction(&cfg, grid);
            V.addPotential(I);
            break;
        default:
            cerr << "Interaction not implemented:: " << interactionType << endl;
            exit(EXIT_FAILURE);
        }
    }

    V.updatePositionBasisElements();
}
//------------------------------------------------------------------------------
void MctdhfApplication::setOneBodyPotentials(SingleParticleOperator &h, const Grid &grid)
{
    Potential* I;
    const Setting& root = cfg.getRoot();
    const Setting &oneBodyPotentials = root["oneBodyPotential"]["potential"];
    int nPotentials = oneBodyPotentials.getLength();

    for(int i=0; i<nPotentials; i++){
        int potential = oneBodyPotentials[i];
        switch (potential) {
        case OB_HARMONIC_OSCILLATOR:
            I = new HarmonicOscillatorOneBody(&cfg, grid);
            h.addPotential(I);
            break;
        case OB_COULOMB_INTERACTION_NUCLEUS:
            I = new CoulombInteractionNucleus(&cfg, grid);
            h.addPotential(I);
            break;
        case OB_ANHARMONIC_DOUBLE_WELL:
            I = new AnharmonicDoubleWell(&cfg, grid);
            h.addPotential(I);
            break;
        case OB_FINITE_HARMONIC_OSCILLATOR:
            I = new FiniteHarmonicOscillator_OB(&cfg, grid);
            h.addPotential(I);
            break;
        case OB_GAUSSIAN_DOUBLE_WELL:
            I = new GaussianDoubleWell(&cfg, grid);
            h.addPotential(I);
            break;
        default:
            cerr << "Potential not implemented:: " << potential << endl;
            exit(EXIT_FAILURE);
        }
    }
}
//------------------------------------------------------------------------------
void MctdhfApplication::setTimeDepOneBodyPotentials(SingleParticleOperator &h, const Grid &grid)
{
    Potential* I;

    const Setting& root = cfg.getRoot();
    const Setting &oneBodyPotentials = root["oneBodyPotential"]["timeDepPotential"];
    int nPotentials = oneBodyPotentials.getLength();

    for(int i=0; i<nPotentials; i++){
        int potential = oneBodyPotentials[i];

        switch (potential) {
        case OB_SIMPLE_LASER:
            I = new simpleLaser(&cfg, grid);
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
        I = new RungeKuttaFehlberg(&cfg);
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
DifferentialOperator* MctdhfApplication::setDifferentialOpertor(const Grid &grid)
{
    // Setting the integrator
    DifferentialOperator *I;
    int differentialOperator = cfg.lookup("spatialDiscretization.differentialOperator");

    switch (differentialOperator) {
    case DO_FINITE_DIFFERENCE_1d:
        I = new FiniteDifference1d(&cfg, grid);
        break;
    case DO_FINITE_DIFFERENCE_FIVE_POINT_1D:
        I = new FiniteDifferenceFivePoint1d(&cfg, grid);
        break;
    case DO_FOURIER_1D:
        I = new Spectral1d(&cfg, grid);
        break;
    case DO_FINITE_DIFFERENCE_2d:
        I = new FiniteDifference2d(&cfg, grid);
        break;
    case DO_FiniteDifferenceFivePoint2d:
        I = new FiniteDifferenceFivePoint2d(&cfg, grid);
        break;
    case DO_FOURIER_2D:
        I = new Fourier2d(&cfg, grid);
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
