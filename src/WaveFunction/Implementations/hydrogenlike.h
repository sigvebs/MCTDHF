#ifndef HYDROGENLIKE_H
#define HYDROGENLIKE_H

#include "../wavefunction.h"

class HydrogenLike: public Wavefunction
{
public:
    HydrogenLike(Config *cfg, vec quantumNumbers);
    virtual double evaluate(const vec &r);
protected:
    int n;
    double laguerrePolynomial(const int n, const double x);
};

#endif // HYDROGENLIKE_H
