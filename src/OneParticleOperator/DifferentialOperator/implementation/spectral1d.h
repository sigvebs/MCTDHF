#ifndef SPECTRAL1D_H
#define SPECTRAL1D_H

// Local includes
#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

// Library includes
//#include <fftw3.h>

class Spectral1d: public DifferentialOperator
{
public:
    Spectral1d(Config *cfg);
    virtual cx_vec secondDerivative(const cx_vec &phi);
};

#endif // SPECTRAL1D_H
