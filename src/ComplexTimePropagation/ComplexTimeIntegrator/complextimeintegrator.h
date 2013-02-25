#ifndef COMPLEXTIMEINTEGRATOR_H
#define COMPLEXTIMEINTEGRATOR_H

// Local includes
#include <src/includes/defines.h>
#include <src/SlaterEquation/slaterequation.h>
#include <src/OrbitalEquation/orbitalequation.h>
#include <src/Interaction/interaction.h>
#include <src/OneParticleOperator/oneparticleoperator.h>

// Library incldues
#include <armadillo>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class ComplexTimeIntegrator
{
public:
    ComplexTimeIntegrator(Config *cfg);
    virtual void stepForward() = 0;
    void setDependencies(SlaterEquation *slater,
                         OrbitalEquation *orbital,
                         Interaction *V,
                         SingleParticleOperator *h);
    cx_vec getCoefficients();
    void setInititalState(cx_vec A, cx_mat C);
protected:
    Config* cfg;
    double dt;
    double t;

    cx_vec A;
    cx_mat C;

    SlaterEquation *slater;
    OrbitalEquation *orbital;

    Interaction *V;
    SingleParticleOperator *h;

    double E;
    vec EVec;
    int step;
    cx_double i;
};

#endif // COMPLEXTIMEINTEGRATOR_H
