#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "src/includes/lib.h"
#include "src/includes/defines.h"

#include <armadillo>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class Wavefunction
{
public:
    Wavefunction(Config *cfg, vec quantumNumbers);
    ~Wavefunction();
    virtual double evaluate(const vec &r) = 0;
protected:
    Config *cfg;
    vec quantumNumbers;
    double coefficient;
    int dim;
};

#endif // WAVEFUNCTION_H
