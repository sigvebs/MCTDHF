#include "finitedifferencefivepoint1d.h"
//------------------------------------------------------------------------------
FiniteDifferenceFivePoint1d::FiniteDifferenceFivePoint1d(Config *cfg):
    DifferentialOperator(cfg)
{
}
//------------------------------------------------------------------------------
cx_vec FiniteDifferenceFivePoint1d::secondDerivative(const cx_vec &phi)
{
    cx_vec diff = zeros<cx_vec>(nGrid);

    for(int j=2; j<nGrid-2; j++){
        diff(j) = - phi(j+2)
                + cx_double(16,0)* phi(j+1)
                - cx_double(30,0)* phi(j)
                + cx_double(16,0)* phi(j-1)
                - phi(j-2);
    }
    // Endpoints
    int j;

    // Start
    j = 1;
    diff(j) = - phi(j+2)
            + cx_double(16,0)* phi(j+1)
            - cx_double(30,0)* phi(j)
            + cx_double(16,0)* phi(j-1);

    j = 0;
    diff(j) = - phi(j+2)
            + cx_double(16,0)* phi(j+1)
            - cx_double(30,0)* phi(j);

    // End
    j = nGrid-2;
    diff(j) = + cx_double(16,0)* phi(j+1)
            - cx_double(30,0)* phi(j)
            + cx_double(16,0)* phi(j-1)
            - phi(j-2);

    j = nGrid-1;
    diff(j) =
            - cx_double(30,0)* phi(j)
            + cx_double(16,0)* phi(j-1)
            - phi(j-2);

    return diff/(12*dx*dx);
}
//------------------------------------------------------------------------------
