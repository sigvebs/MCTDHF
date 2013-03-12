#ifndef POTENTIAL_H
#define POTENTIAL_H

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

class Potential
{
public:
    Potential(Config* cfg);
    virtual cx_vec evaluate(const cx_vec &psi, double t=0) = 0;
protected:
    Config* cfg;

    vec potential;
    int nGrid;
};

#endif // POTENTIAL_H
