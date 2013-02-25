#include "complextimeintegrator.h"

//------------------------------------------------------------------------------
ComplexTimeIntegrator::ComplexTimeIntegrator(Config *cfg)
{
    try {
        dt = cfg->lookup("ComplexTimeIntegration.dt");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ComplexTimeIntegrator::ComplexTimeIntegrator(Config *cfg)::"
             << "Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
    E = 0;
    step = 0;
    t = 0;
    int N = cfg->lookup("ComplexTimeIntegration.N");
    EVec = zeros(N);
#if DEBUG
    cout << "ComplexTimeIntegrator::ComplexTimeIntegrator(Config *cfg)" << endl;
    cout << "dt = " << dt << endl;
#endif
}
//------------------------------------------------------------------------------
void ComplexTimeIntegrator::setDependencies(SlaterEquation *slater,
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
cx_vec ComplexTimeIntegrator::getCoefficients()
{
    return A;
}
//------------------------------------------------------------------------------
void ComplexTimeIntegrator::setInititalState(cx_vec A, cx_mat C)
{
    this->A = A;
    this->C = C;
}
//------------------------------------------------------------------------------
