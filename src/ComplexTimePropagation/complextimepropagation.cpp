#include "complextimepropagation.h"

//------------------------------------------------------------------------------
ComplexTimePropagation::ComplexTimePropagation(Config *cfg):
    cfg(cfg)
{
    try{
        dt = cfg->lookup("ComplexTimeIntegration.dt");
        cfg->lookupValue("systemSettings.filePath", filePath);
        N = cfg->lookup("ComplexTimeIntegration.N");
        saveToFileInterval = cfg->lookup("systemSettings.saveToFileInterval");
        printProgress  = cfg->lookup("systemSettings.printProgress");
        saveEveryTimeStep = cfg->lookup("systemSettings.saveEveryTimeStep");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg)::Error reading parameter from config object." << endl;
        exit(EXIT_FAILURE);
    }
    step = 0;
    t = 0;
    E = zeros(N/saveToFileInterval+1);
    dE = zeros(N/saveToFileInterval+1);
    K = vec(1);
    Eprev = 99;
    filenameOrbitals = filePath + "C.mat";
    filenameSlaterDet = filePath + "A.mat";
    filenameEnergy = filePath + "E.mat";
    filenameDeltaE = filePath + "dE.mat";
    filenameSvdRho = filePath + "svdRho.mat";
    filenameRho = filePath + "rho.mat";
    filenameCorrelation = filePath + "K.mat";
#ifdef DEBUG
    cout << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg)::" << endl
         << "dt \t= " << dt << endl
         << "N \t= " << N << endl;
#endif
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::doComplexTimePropagation()
{
    cout.precision(16);
    int counter = 0;
    bool accepted;

    cout << "Starting imaginary propagation:" << endl;
    for(step=0; step < N; step++){
        accepted = this->stepForward();

        // Saving C and A to disk
        if((step % saveToFileInterval == 0 || step == N-1) && accepted){

            // Updating the one-body- and interaction-elements
            V->computeNewElements(C);
            h->computeNewElements(C);

            // Collecting data
            E(counter) = slater->getEnergy(A) ;
            dE(counter) = E(counter) - Eprev;
            rho = &orbital->reCalculateRho1(A);
            K = orbital->getCorrelation();
            svdRho = orbital->getSvdRho1();

            saveProgress(counter);
            printProgressToScreen(counter);

            Eprev = E(counter);
            counter++;
        }
    }
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::setDependencies(SlaterEquation *slater,
                                            OrbitalEquation *orbital,
                                            Interaction *V,
                                            SingleParticleOperator *h)
{
    this->slater = slater;
    this->orbital = orbital;
    this->V = V;
    this->h = h;
}
//------------------------------------------------------------------------------
cx_vec ComplexTimePropagation::getCoefficients()
{
    return A;
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::setInititalState(cx_vec &A, cx_mat &C)
{
    this->A = A;
    this->C = C;
    nOrbitals = C.n_cols;
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::renormalize(cx_mat &D)
{
    // Re-normalization of C using SVD
    cx_mat X;
    vec s;
    cx_mat Y;
    svd_econ(X, s, Y, D);
    D = X*Y.t();
}
//------------------------------------------------------------------------------
cx_mat ComplexTimePropagation::getCurrentC()
{
    return C;
}
//------------------------------------------------------------------------------
cx_vec ComplexTimePropagation::getCurrentA()
{
    return A;
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::printProgressToScreen(uint counter)
{
    if( printProgress ){
        cout << "---------------------------------------------------------------\n";
        cout << "step = " << step << endl
             << "E = " << E(counter) << endl
             << "dE = " << dE(counter) << endl
             << "k = " << K(0) << endl
             << "dt = " << dt << endl
             << "svdRho1 = ";
        svdRho.raw_print();
        cout << endl;
    }
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::saveProgress(uint counter)
{
    h->saveOperators();
    C.save(filenameOrbitals);
    A.save(filenameSlaterDet, arma_ascii);
    E.save(filenameEnergy);
    dE.save(filenameDeltaE);
    K.save(filenameCorrelation);
    svdRho.save(filenameSvdRho);
    (*rho).save(filenameRho);

    if(saveEveryTimeStep){
        stringstream fileName;
        fileName << filePath << "cx_C" << counter << ".mat";
        C.save(fileName.str());
        fileName.str("");
        fileName << filePath << "cx_A" << counter << ".mat";
        A.save(fileName.str());
        fileName.str("");
        fileName << filePath << "cx_rho" << counter << ".mat";
         (*rho).save(fileName.str());
    }
}
//------------------------------------------------------------------------------
