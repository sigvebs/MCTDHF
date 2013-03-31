#include "finitedifference2d.h"

//------------------------------------------------------------------------------
FiniteDifference2d::FiniteDifference2d(Config *cfg, const Grid &grid):
    DifferentialOperator(cfg, grid)
{
    nGridX = grid.nGridX;
    nGridY = grid.nGridY;

    diff = zeros<cx_vec>(nGrid);
}
//------------------------------------------------------------------------------
cx_vec FiniteDifference2d::secondDerivative(const cx_vec &phi)
{
    diff.zeros();

    cx_double diffX;
    cx_double diffY;

    for(int i=0; i<nGrid; i++){

        // Center point
        diffX = -2*phi(i);
        diffY = diffX;


        // X -----------------------
        if(i + nGridX < nGrid)
            diffX += phi(i + nGridX);

        if(i - nGridX > 0)
            diffX += phi(i - nGridX);

        // Y -----------------------
        if(i + 1 < nGrid)
            diffY += phi(i + 1);

        if(i - 1 > 0)
            diffY += phi(i - 1);

        diff(i) = diffX/dxdx + diffY/dydy;

    }

    return diff;
}
//------------------------------------------------------------------------------
