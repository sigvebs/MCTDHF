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
    cx_mat computeRightHandSide(const cx_mat &C, const cx_vec &A);
protected:
    void computeProjector(const cx_mat &C);
    void computeUMatrix(const cx_mat &C);
    void computeOneParticleReducedDensity();
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
    int dim;
    double dx;

    cx_mat rho;
    unordered_map<int, cx_double> rho2;
    const cx_vec* A;
    const cx_mat* hC;
    cx_mat U;
    cx_mat P;
    cx_mat I ;

    cx_mat rightHandSide;

    // REMOVE
    int counter;
};

#endif // ORBITALEQUATION_H
