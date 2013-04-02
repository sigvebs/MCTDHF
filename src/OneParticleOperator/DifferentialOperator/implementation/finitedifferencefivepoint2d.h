#ifndef FINITEDIFFERENCEFIVEPOINT2D_H
#define FINITEDIFFERENCEFIVEPOINT2D_H

// Local includes
#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

class FiniteDifferenceFivePoint2d: public DifferentialOperator
{
public:
    FiniteDifferenceFivePoint2d(Config* cfg, const Grid &grid);
    virtual cx_vec secondDerivative(const cx_vec &phi);
protected:
    cx_vec diff;
    int nGridX;
    int nGridY;
};

#endif // FINITEDIFFERENCEFIVEPOINT2D_H
