#ifndef SPECTRAL1D_H
#define SPECTRAL1D_H

// Local includes
#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

// Library includes
#include <fftw3.h>

class Spectral1d: public DifferentialOperator
{
public:
    Spectral1d(Config *cfg, const vec &x);
    ~Spectral1d();
    virtual cx_vec secondDerivative(const cx_vec &phi);
protected:
    vec k;
    cx_vec diff;
    fftw_plan forward;
    fftw_plan backward;
    double scaling;
};

#endif // SPECTRAL1D_H
