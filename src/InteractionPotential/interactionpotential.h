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
    InteractionPotential(Config* cfg);
    virtual mat computeInteractionSpace() = 0;
protected:
    Config* cfg;

    mat interactionSpace;
    vec x;
    int nGrid;
};

#endif // INTERACTIONPOTENTIAL_H
