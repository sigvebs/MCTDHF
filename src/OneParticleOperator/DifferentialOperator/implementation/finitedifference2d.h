#ifndef FINITEDIFFERENCE2D_H
#define FINITEDIFFERENCE2D_H

// Local includes
#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

class FiniteDifference2d: public DifferentialOperator
{
public:
    FiniteDifference2d(Config* cfg, const Grid &grid);
    virtual cx_vec secondDerivative(const cx_vec &phi);
protected:
    cx_vec diff;
    int nGridX;
    int nGridY;
};

#endif // FINITEDIFFERENCE2D_H
