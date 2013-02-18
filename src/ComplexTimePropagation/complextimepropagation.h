#ifndef COMPLEXTIMEPROPAGATION_H
#define COMPLEXTIMEPROPAGATION_H

// Local includes
#include <src/includes/defines.h>
#include <src/OrbitalEquation/orbitalequation.h>
#include <src/SlaterEquation/slaterequation.h>
#include <src/ComplexTimePropagation/ComplexTimeIntegrator/complextimeintegrator.h>

// Libary incldues
#include <armadillo>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;
using namespace arma;

class ComplexTimePropagation
{
public:
    ComplexTimePropagation(Config *cfg, ComplexTimeIntegrator* complexTimeIntegrator);
    void doComplexTimePropagation();
protected:
    Config *cfg;
    ComplexTimeIntegrator* complexTimeIntegrator;

    double dt;
    int N;
};

#endif // COMPLEXTIMEPROPAGATION_H
