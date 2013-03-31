#include "finitedifferencefivepoint2d.h"

//------------------------------------------------------------------------------
FiniteDifferenceFivePoint2d::FiniteDifferenceFivePoint2d(Config *cfg, const Grid &grid):
    DifferentialOperator(cfg, grid)
{
    nGridX = grid.nGridX;
    nGridY = grid.nGridY;

    diff = zeros<cx_vec>(nGrid);
}
//------------------------------------------------------------------------------
cx_vec FiniteDifferenceFivePoint2d::secondDerivative(const cx_vec &phi)
{
    diff.zeros();

    cx_double diffX;
    cx_double diffY;

    for(int i=0; i<nGrid; i++){

        // Center point
        diffX = -(cx_double)30*phi(i);
        diffY = diffX;

        // X -----------------------
        if(i + nGridX < nGrid)
            diffX += (cx_double)16*phi(i + nGridX);

        if(i - nGridX > 0)
            diffX += (cx_double)16*phi(i - nGridX);

        if(i + 2*nGridX < nGrid)
            diffX -= phi(i + 2*nGridX);

        if(i - 2*nGridX > 0)
            diffX -= phi(i - 2*nGridX);

        // Y -----------------------
        if(i + 1 < nGrid)
            diffY += (cx_double)16*phi(i + 1);

        if(i - 1 > 0)
            diffY += (cx_double)16*phi(i - 1);

        if(i + 2 < nGrid)
            diffY -= phi(i + 2);

        if(i - 2 > 0)
            diffY -= phi(i - 2);

        diff(i) = diffX/(12*dxdx) + diffY/(12*dydy);
    }

    return diff;
}
//------------------------------------------------------------------------------
