#ifndef SLATEREQUATION_H
#define SLATEREQUATION_H

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

using namespace libconfig;
using namespace std;
using namespace arma;

class SlaterEquation
{
public:
    SlaterEquation(Config *cfg,
                   vector<vec> orbitals,
                   vector<bitset<BITS> > slaterDeterminants,
                   Interaction *interaction,
                   SingleParticleOperator *oneParticleOperator);
    ~SlaterEquation();
    cx_vec computeRightHandSide(const cx_vec &A);
    cx_vec computeRightHandSideComplexTime(const cx_vec &A);
    double getEnergy(const cx_vec &A);
protected:
    void computeHamiltonianMatrix();
    cx_double secondQuantizationOneBodyOperator(const int p, const int q,
                                                bitset<BITS> state1,
                                                const bitset<BITS> &state2);
    cx_double secondQuantizationTwoBodyOperator(const int p,const int q,
                                                const int r, const int s,
                                                bitset<BITS> state1,
                                                const bitset<BITS> &state2);
    // Class variables
    Config* cfg;

    vector<bitset<BITS> > slaterDeterminants;
    vector<vec> orbitals;

    cx_mat H;
    const cx_mat* h;
    int nSlaterDeterminants;
    int nSpatialOrbitals;
    int nOrbitals;
    int dim;

    Interaction* interaction;
    SingleParticleOperator* oneParticleOperator;
};

#endif // SLATEREQUATION_H
