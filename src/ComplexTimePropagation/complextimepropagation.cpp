#include "complextimepropagation.h"

//------------------------------------------------------------------------------
ComplexTimePropagation::ComplexTimePropagation(Config *cfg):
    cfg(cfg)
{
    string filePath;
    try{
        dt = cfg->lookup("ComplexTimeIntegration.dt");
        cfg->lookupValue("systemSettings.filePath", filePath);
        N = cfg->lookup("ComplexTimeIntegration.N");
        cfg->lookupValue("systemSettings.filePath", filePath);
        saveToFileInterval = cfg->lookup("systemSettings.saveToFileInterval");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg)::Error reading parameter from config object." << endl;
        exit(EXIT_FAILURE);
    }
    step = 0;
    t = 0;
    E = zeros(N/saveToFileInterval+1);
    filenameOrbitals = filePath + "C.mat";
    filenameSlaterDet = filePath + "A.mat";
    filenameEnergy = filePath + "E.mat";
#ifdef DEBUG
    cout << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg)::" << endl
         << "dt \t= " << dt << endl
         << "N \t= " << N << endl;
#endif
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::doComplexTimePropagation()
{
    double correlation;
    cout.precision(16);
    int counter = 0;
    bool accepted;
    cout << "Starting imaginary propagation:" << endl;
    for(step=0; step < N; step++){

        accepted = this->stepForward();

        // Saving C and A to disk
        if((step % saveToFileInterval == 0 || step == N-1) && accepted){
            E(counter) = slater->getEnergy(A) ;
            correlation = orbital->getCorrelation(A);
            h->saveOperators();
            cout << "---------------------------------------------------------------\n";
            C.save(filenameOrbitals, arma_ascii);
            A.save(filenameSlaterDet, arma_ascii);
            E.save(filenameEnergy, arma_ascii);
            cout << "step = " << step << endl
                 << "E = " << E(counter) << endl
                 << "k = " << correlation << endl
                 << "dt = " << dt << endl;
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
void ComplexTimePropagation::setInititalState(cx_vec A, cx_mat C)
{
    this->A = A;
    this->C = C;
}
//------------------------------------------------------------------------------
