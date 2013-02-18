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
                    vector<vec> orbitals,
                    vector<bitset<BITS> > slaterDeterminants,
                    Interaction *V,
                    SingleParticleOperator *h);
    void setInititalState(cx_mat C);
    void setSlaterCoefficients(const cx_vec &A);
    const cx_mat &getCoefficientMatrix() const;
    cx_mat computeRightHandSide(cx_mat C, cx_vec A);
protected:
    void computeOneBodyMatrix();
    void computeProjector();
    void computeUMatrix();
    void computeOneParticleReducedDensity();
    void computeTwoParticleReducedDensity();
    cx_double reducedTwoParticleOperator(int p, int q, int s, int r);

    // Class variables
    Config* cfg;

    Interaction *V;
    SingleParticleOperator* h;

    vector<vec> orbitals;
    vector<bitset<BITS> > slaterDeterminants;
    unordered_map<int, cx_double> rho2;

    int nOrbitals;
    int nSpatialOrbitals;
    int nGrid;
    int nSlaterDeterminants;
    int nParticles;
    int dim;

    cx_mat rho;
    cx_mat C;
    cx_vec A;
    cx_mat U;
    cx_mat P;
    cx_mat hC;
    cx_mat I ;
};

#endif // ORBITALEQUATION_H
