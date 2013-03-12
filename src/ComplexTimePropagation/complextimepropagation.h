#ifndef COMPLEXTIMEPROPAGATION_H
#define COMPLEXTIMEPROPAGATION_H

// Local includes
#include <src/includes/defines.h>
#include <src/OrbitalEquation/orbitalequation.h>
#include <src/SlaterEquation/slaterequation.h>

// Libary incldues
#include <armadillo>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;
using namespace arma;

class ComplexTimePropagation
{
public:
    ComplexTimePropagation(Config *cfg);
    void doComplexTimePropagation();

    virtual bool stepForward() = 0;
    void setDependencies(SlaterEquation *slater,
                         OrbitalEquation *orbital,
                         Interaction *V,
                         SingleParticleOperator *h);
    cx_vec getCoefficients();
    void setInititalState(cx_vec A, cx_mat C);
protected:
    Config *cfg;
    double dt;
    double t;

    cx_vec A;
    cx_mat C;

    SlaterEquation *slater;
    OrbitalEquation *orbital;

    Interaction *V;
    SingleParticleOperator *h;

    vec E;
    int step;
    cx_double i;
    int N;

    // Filenames
    int saveToFileInterval;
    string filenameOrbitals;
    string filenameSlaterDet;
    string filenameEnergy;
};

#endif // COMPLEXTIMEPROPAGATION_H
