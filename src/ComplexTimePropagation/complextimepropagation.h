#ifndef COMPLEXTIMEPROPAGATION_H
#define COMPLEXTIMEPROPAGATION_H

// Local includes
#include <src/includes/defines.h>
#include <src/OrbitalEquation/orbitalequation.h>
#include <src/SlaterEquation/slaterequation.h>

// Libary incldues
#include <armadillo>
#include <libconfig.h++>

#ifdef USE_MPI
#include <mpi.h>
#endif

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
    void setInititalState(cx_vec &A, cx_mat &C);
    void renormalize(cx_mat &D);
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

    bool printProgress;
    bool saveEveryTimeStep;

    int myRank, nNodes;
    bool isMaster;
};

#endif // COMPLEXTIMEPROPAGATION_H
