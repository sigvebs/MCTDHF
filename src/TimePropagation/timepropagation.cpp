#include "timepropagation.h"

//------------------------------------------------------------------------------
TimePropagation::TimePropagation(Config *cfg):
    cfg(cfg)
{
    string filePath;
    try{
        dt = cfg->lookup("timeIntegration.dt");
        cfg->lookupValue("systemSettings.filePath", filePath);
        N = cfg->lookup("timeIntegration.N");
        cfg->lookupValue("systemSettings.filePath", filePath);
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg)::Error reading parameter from config object." << endl;
        exit(EXIT_FAILURE);
    }
    step = 0;
    t = 0;
    overlap = zeros(N);
    time = zeros(N);
    filenameOverlap = filePath + "overlap.mat";
    filenameT = filePath + "time.mat";
#ifdef DEBUG
    cout << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg)::" << endl
         << "dt \t= " << dt << endl
         << "N \t= " << N << endl;
#endif
}
//------------------------------------------------------------------------------
void TimePropagation::doComplexTimePropagation()
{
    cout << "\n" << endl;
    bool accepted;

    for(step=0; step < N; step++){
        accepted = this->stepForward();

        // Saving C and A to disk
        if(accepted){
            time(step) = t;
            overlap(step) = pow(fabs(cdot(A0, A)),2);
//            cout << "-------------------------------------------------------\n";
//            cout << "step = " << step << endl
//                 << "E = " << slater->getEnergy(A) << endl
//                 << "t = " << t << endl
//                 << "|<A0|A>|^2 = " << overlap(step) << endl
//                 << "dt = " << dt << endl;
            t += dt;
        }
    }

    // Saving results
    time.save(filenameT, arma_ascii);
    overlap.save(filenameOverlap, arma_ascii);
}
//------------------------------------------------------------------------------
void TimePropagation::setDependencies(SlaterEquation *slater,
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
void TimePropagation::setInititalState(cx_vec A, cx_mat C)
{
    this->A = A;
    A0 = A;
    this->C = C;
    C0 = C;
}
//------------------------------------------------------------------------------
double TimePropagation::computeOverlap()
{
    return 0;
}
//------------------------------------------------------------------------------

