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
    virtual mat computeInteractionSpace() = 0;
protected:
    Config* cfg;
    mat interactionSpace;
    vec ri;
    vec rj;

    const Grid &grid;
    int nGrid;
    int dim;
};

#endif // INTERACTIONPOTENTIAL_H
