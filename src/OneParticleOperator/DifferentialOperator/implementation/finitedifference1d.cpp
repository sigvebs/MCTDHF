#include "finitedifference1d.h"

//------------------------------------------------------------------------------
FiniteDifference1d::FiniteDifference1d(Config *cfg, const Grid &grid):
    DifferentialOperator(cfg, grid)
{
}
//------------------------------------------------------------------------------
cx_vec FiniteDifference1d::secondDerivative(const cx_vec &phi)
{
    cx_vec diff = zeros<cx_vec>(nGrid);

    for(int j=1; j<nGrid-1; j++){
        diff(j) = phi(j-1) - cx_double(2,2)*phi(j) + phi(j+1);
    }

    // Endpoints
    diff(0) = - cx_double(2,2)*phi(0) + phi(1);
    diff(nGrid-1) = phi(nGrid-2) - cx_double(2,2)*phi(nGrid-1);

    return diff/(dxdx);
}
//------------------------------------------------------------------------------
