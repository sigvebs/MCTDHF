#include "finitedifference1d.h"
//------------------------------------------------------------------------------
FiniteDifference1d::FiniteDifference1d(Config *cfg):
    DifferentialOperator(cfg)
{
}
//------------------------------------------------------------------------------
cx_vec FiniteDifference1d::secondDerivative(const cx_vec &phi)
{
    cx_vec diff = zeros<cx_vec>(nGrid);

    for(int j=1; j<nGrid-1; j++){
        diff(j) = phi(j-1) - cx_double(2,0)*phi(j) + phi(j+1);
    }
    // Endpoints
    //    diff(0) = phi(nGrid-1) - cx_double(2,0)*phi(0) + phi(1);
    //    diff(nGrid-1) = phi(nGrid-2) - cx_double(2,0)*phi(nGrid-1) + phi(0);

    diff(0) = - cx_double(2,0)*phi(0) + phi(1);
    diff(nGrid-1) = phi(nGrid-2) - cx_double(2,0)*phi(nGrid-1);

    //    diff(0) = 0;
    //    diff(nGrid-1) = 0;

    return diff/(dx*dx);
}
//------------------------------------------------------------------------------
