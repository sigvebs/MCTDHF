#ifndef HARMONICOSCILLATOR_H
#define HARMONICOSCILLATOR_H

#include "../wavefunction.h"

class HarmonicOscillator: public Wavefunction
{
public:
    HarmonicOscillator(Config *cfg, vec quantumNumbers);
    virtual mat evaluate(const mat &x) = 0;
protected:
    double hermitePolynomial(const int degree, const double x);
    double w;
    double sqrtW;
};

#endif // HARMONICOSCILLATOR_H
