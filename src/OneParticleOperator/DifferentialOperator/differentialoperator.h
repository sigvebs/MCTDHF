#ifndef DIFFERENTIALOPERATOR_H
#define DIFFERENTIALOPERATOR_H

// Local includes
#include <src/includes/defines.h>
#include <src/includes/lib.h>

// Library includes
#include <armadillo>
#include <vector>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class DifferentialOperator
{
public:
    DifferentialOperator(Config* cfg, const vec &x);
    virtual cx_vec secondDerivative(const cx_vec &phi) = 0;
protected:
    Config* cfg;

    int nGrid;
    double dx;
    double dxdx;
    vec x;
};

#endif // DIFFERENTIALOPERATOR_H
