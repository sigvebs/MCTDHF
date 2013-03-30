#ifndef HYDROGENLIKE_H
#define HYDROGENLIKE_H

#include "../wavefunction.h"

class HydrogenLike: public Wavefunction
{
public:
    HydrogenLike(Config *cfg, vec quantumNumbers);
    virtual double evaluate(double x);
    virtual double evaluate(double x, double y);
protected:
    int n;
    double laguerrePolynomial(const int n, const double x);
};

#endif // HYDROGENLIKE_H
