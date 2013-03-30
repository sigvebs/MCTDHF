#ifndef HARMONICOSCILLATOR_H
#define HARMONICOSCILLATOR_H

#include "../wavefunction.h"

class HarmonicOscillator: public Wavefunction
{
public:
    HarmonicOscillator(Config *cfg, vec quantumNumbers);
    virtual double evaluate(double x);
    virtual double evaluate(double x, double y);
protected:
    double hermitePolynomial(const int degree, const double x);
    double w;
    double sqrtW;
};

#endif // HARMONICOSCILLATOR_H
