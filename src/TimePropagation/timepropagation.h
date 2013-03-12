#ifndef TIMEPROPAGATION_H
#define TIMEPROPAGATION_H


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

class TimePropagation
{
public:
    TimePropagation(Config *cfg);
    void doComplexTimePropagation();

    virtual bool stepForward() = 0;
    void setDependencies(SlaterEquation *slater,
                         OrbitalEquation *orbital,
                         Interaction *V,
                         SingleParticleOperator *h);
    void setInititalState(cx_vec A, cx_mat C);
protected:
    double computeOverlap();
    Config *cfg;
    double dt;
    double t;
    vec time;

    cx_vec A;
    cx_vec A0;
    cx_mat C;
    cx_mat C0;

    SlaterEquation *slater;
    OrbitalEquation *orbital;

    Interaction *V;
    SingleParticleOperator *h;

    vec overlap;
    int step;
    cx_double i;
    int N;

    // Filenames
    string filenameOverlap;
    string filenameT;
};

#endif // TIMEPROPAGATION_H
