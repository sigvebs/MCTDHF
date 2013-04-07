#ifndef FOURIER2D_H
#define FOURIER2D_H

// Local includes
#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

// Library includes
#include <fftw3.h>

class Fourier2d: public DifferentialOperator
{
public:
    Fourier2d(Config *cfg, const Grid &grid);
    ~Fourier2d();
    virtual cx_vec secondDerivative(const cx_vec &phi);
protected:
    vec k;
    cx_vec diff;
    fftw_plan forward;
    fftw_plan backward;

    int nGridX;
    int nGridY;
};

#endif // FOURIER2D_H
