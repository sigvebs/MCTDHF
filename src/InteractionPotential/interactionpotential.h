#ifndef INTERACTIONPOTENTIAL_H
#define INTERACTIONPOTENTIAL_H
// Library includes
#include <armadillo>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class InteractionPotential
{
public:
    InteractionPotential(Config* cfg, const vec &x);
    virtual mat computeInteractionSpace() = 0;
protected:
    Config* cfg;

    mat interactionSpace;
    const vec &x;
    int nGrid;
};

#endif // INTERACTIONPOTENTIAL_H
