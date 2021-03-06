#ifndef FINITEDIFFERENCE1D_H
#define FINITEDIFFERENCE1D_H

// Local includes
#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

class FiniteDifference1d: public DifferentialOperator
{
public:
    FiniteDifference1d(Config* cfg, const Grid &grid);
    virtual cx_vec secondDerivative(const cx_vec &phi);
};

#endif // FINITEDIFFERENCE1D_H
