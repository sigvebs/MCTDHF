#ifndef FINITEDIFFERENCE1D_H
#define FINITEDIFFERENCE1D_H

#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

class FiniteDifference1d: public DifferentialOperator
{
public:
    FiniteDifference1d(Config* cfg);
    virtual cx_vec secondDerivative(const cx_vec &phi);
};

#endif // FINITEDIFFERENCE1D_H
