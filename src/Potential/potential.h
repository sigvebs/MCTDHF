#ifndef POTENTIAL_H
#define POTENTIAL_H

// Library includes
#include <armadillo>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class Potential
{
public:
    Potential(Config* cfg);
    virtual mat computeInteractionSpace() = 0;
protected:
    Config* cfg;

    mat interactionSpace;
    vec x;
    int nGrid;
};

#endif // POTENTIAL_H
