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
    t = 0;
#if DEBUG
    cout << "ComplexTimeIntegrator::ComplexTimeIntegrator(Config *cfg)" << endl;
    cout << "dt = " << dt << endl;
#endif
}

void ComplexTimeIntegrator::setDependencies(SlaterEquation *slater,
                                            OrbitalEquation *orbital,
                                            Interaction *V,
                                            SingleParticleOperator *h)
{
    this->slater = slater;
    this->orbital = orbital;
    this->V = V;
    this->h = h;
    A = slater->getCoefficientVector();
    C = orbital->getCoefficientMatrix();
}
//------------------------------------------------------------------------------
cx_vec ComplexTimeIntegrator::getCoefficients()
{
    return A;
}
//------------------------------------------------------------------------------
