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
        cerr << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg, SlaterEquation *slaterEquation, OrbitalEquation *orbEq)::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
    step = 0;
    t = 0;
    E = zeros(N);
    filenameOrbitals = filePath + "C.mat";
    fileNameSlaterDet = filePath + "A.mat";
    fileNameEnergy = filePath + "E.mat";
#ifdef DEBUG
    cout << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg, SlaterEquation *slaterEquation, OrbitalEquation *orbEq)" << endl
         << "dt \t= " << dt << endl
         << "N \t= " << N << endl;
#endif
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::doComplexTimePropagation()
{
    for(step=0; step<N; step++){

        this->stepForward();

        // Saving C and A to disk
        if(step % saveToFileInterval == 0 || step == N-1){
            E(step) = slater->getEnergy(A) ;
            cout << "---------------------------------------------------------------\n";
            C.save(filenameOrbitals, arma_ascii);
            A.save(fileNameSlaterDet, arma_ascii);
            E.save(fileNameEnergy, arma_ascii);
            cout.precision(16);
            cout << "step = " << step << endl
                 << "E = " << E(step) << endl;
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
