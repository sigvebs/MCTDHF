#ifndef FINITEDIFFERENCEFIVEPOINT1D_H
#define FINITEDIFFERENCEFIVEPOINT1D_H

#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

class FiniteDifferenceFivePoint1d: public DifferentialOperator
{
public:
    FiniteDifferenceFivePoint1d(Config* cfg);
    virtual cx_vec secondDerivative(const cx_vec &phi);
};

#endif // FINITEDIFFERENCEFIVEPOINT1D_H
