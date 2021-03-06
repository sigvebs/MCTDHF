#ifndef INTERACTIONPOTENTIAL_H
#define INTERACTIONPOTENTIAL_H

//Local includes
#include <src/Grid/grid.h>

// Library includes
#include <armadillo>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class InteractionPotential
{
public:
    InteractionPotential(Config* cfg, const Grid &grid);
    virtual double evaluate(uint i, uint j) = 0;
protected:
    const Grid &grid;
    int nGrid;
    int dim;

    Config* cfg;
};

#endif // INTERACTIONPOTENTIAL_H
