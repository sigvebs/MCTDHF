#ifndef DIFFERENTIALOPERATOR_H
#define DIFFERENTIALOPERATOR_H

// Local includes
#include <src/includes/defines.h>
#include <src/includes/lib.h>
#include <src/Grid/grid.h>

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
    DifferentialOperator(Config* cfg, const Grid &grid);
    virtual cx_vec secondDerivative(const cx_vec &phi) = 0;
protected:
    double dx;
    double dxdx;
    const Grid &grid;
    int nGrid;

    Config* cfg;
};

#endif // DIFFERENTIALOPERATOR_H
