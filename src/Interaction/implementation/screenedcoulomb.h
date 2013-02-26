#ifndef SCREENEDCOULOMB_H
#define SCREENEDCOULOMB_H

// Library includes
#include <armadillo>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class ScreenedCoulomb
{
public:
    ScreenedCoulomb(Config* cfg);
    mat computeInteractionSpace();
protected:
    Config* cfg;

    mat interactionSpace;
    vec x;
    int nGrid;
    double aa;
};

#endif // SCREENEDCOULOMB_H
