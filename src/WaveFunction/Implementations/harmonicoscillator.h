#ifndef HARMONICOSCILLATOR_H
#define HARMONICOSCILLATOR_H

#include "../wavefunction.h"

class HarmonicOscillator: public Wavefunction
{
public:
    HarmonicOscillator(Config *cfg, vec quantumNumbers);
    virtual mat evaluate(const mat &x) = 0;
    virtual double getEnergy() = 0;
protected:
    double hermitePolynomial(const int degree, const double x);
    int dim;
    double w;
    double sqrtW;
    double aa;
};

#endif // HARMONICOSCILLATOR_H
