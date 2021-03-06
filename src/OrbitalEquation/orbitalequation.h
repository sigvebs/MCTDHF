#ifndef ORBITALEQUATION_H
#define ORBITALEQUATION_H

// Local includes
#include <src/includes/defines.h>
#include <src/includes/lib.h>
#include <src/includes/binaryoperations.h>

#include <src/Interaction/interaction.h>
#include <src/OneParticleOperator/oneparticleoperator.h>

// Libraries
#include <libconfig.h++>
#include <armadillo>
#include <bitset>
#include <unordered_map>

using namespace libconfig;
using namespace std;
using namespace arma;

class OrbitalEquation
{
public:
    OrbitalEquation(Config* cfg,
                    vector<bitset<BITS> > slaterDeterminants,
                    Interaction *V,
                    SingleParticleOperator *h);
    const cx_mat &computeRightHandSide(const cx_mat &C, const cx_vec &A);
    double getCorrelation();
    const cx_mat &reCalculateRho1(const cx_vec &A);
    vec getSvdRho1();
protected:
    void computeProjector(const cx_mat &C);
    cx_double findRho2(int p,int q, int r, int s);
    void computeUMatrix(const cx_mat &C);
    void computeOneParticleReducedDensity();
    void computeOneParticleReducedDensityWithSpin();
    void computeTwoParticleReducedDensity();
    cx_double reducedOneParticleOperator(const int i, const int j);
    cx_double reducedTwoParticleOperator(const int p, const int q,
                                         const int r, const int s);

    // Class variables
    Config* cfg;

    Interaction *V;
    SingleParticleOperator* h;

    vector<bitset<BITS> > slaterDeterminants;

    int nOrbitals;
    int nGrid;
    int nSlaterDeterminants;
    int nParticles;

    cx_mat invRho;
    unordered_map<int, cx_double> rho2;
    const cx_vec* A;
    const cx_mat* hC;
    cx_mat U;
    cx_mat Q;

    cx_mat rightHandSide;
    cx_mat rho1; // TMP

    // MPI
    imat allRH;
    int myRank, nNodes;
    vector<pair<int,int> > myRij;
    ivec sizeRij;
};

#endif // ORBITALEQUATION_H
