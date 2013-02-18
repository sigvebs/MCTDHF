#include "complextimepropagation.h"

//------------------------------------------------------------------------------
ComplexTimePropagation::ComplexTimePropagation(Config *cfg, ComplexTimeIntegrator *complexTimeIntegrator):
    cfg(cfg),
    complexTimeIntegrator(complexTimeIntegrator)
{
    try{
        dt = cfg->lookup("ComplexTimeIntegration.dt");
        N = cfg->lookup("ComplexTimeIntegration.N");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg, SlaterEquation *slaterEquation, OrbitalEquation *orbEq)::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
#ifdef DEBUG
    cout << "ComplexTimePropagation::ComplexTimePropagation(Config *cfg, SlaterEquation *slaterEquation, OrbitalEquation *orbEq)" << endl
         << "dt \t= " << dt << endl
         << "N \t= " << N << endl;
#endif
}
//------------------------------------------------------------------------------
void ComplexTimePropagation::doComplexTimePropagation()
{
    for(int i=0; i<N; i++){
        complexTimeIntegrator->stepForward();
    }
}
//------------------------------------------------------------------------------
