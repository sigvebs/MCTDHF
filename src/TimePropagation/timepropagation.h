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
    void doTimePropagation();

    virtual bool stepForward() = 0;
    void setDependencies(SlaterEquation *slater,
                         OrbitalEquation *orbital,
                         Interaction *V,
                         SingleParticleOperator *h);
    void setInititalState(cx_vec &A, cx_mat &C);
    cx_mat getCurrentC();
    cx_vec getCurrentA();
protected:
    void printProgressToScreen(uint counter);
    void saveProgress(uint counter);
    Config *cfg;
    double dt;
    double t;
    int nOrbitals;

    cx_vec A;
    cx_mat C;

    SlaterEquation *slater;
    OrbitalEquation *orbital;

    Interaction *V;
    SingleParticleOperator *h;

    vec E;
    vec dE;
    vec K;
    double Eprev;
    int step;
    cx_double i;
    int N;
    vec svdRho;
    const cx_mat *rho;
    vec time;

    // Filenames
    int saveToFileInterval;
    string filePath;
    string filenameOrbitals;
    string filenameSlaterDet;
    string filenameEnergy;
    string filenameDeltaE;
    string filenameSvdRho;
    string filenameRho;
    string filenameCorrelation;
    string filenameT;

    bool printProgress;
    bool saveEveryTimeStep;
};

#endif // TIMEPROPAGATION_H
