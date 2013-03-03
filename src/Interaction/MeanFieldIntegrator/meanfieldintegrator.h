#ifndef MEANFIELDINTEGRATOR_H
#define MEANFIELDINTEGRATOR_H

// Local includes
#include <src/includes/defines.h>
#include <src/includes/lib.h>
#include <src/InteractionPotential/interactionpotential.h>

// Library includes
#include <armadillo>
#include <vector>
#include <unordered_map>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class MeanFieldIntegrator
{
public:
    MeanFieldIntegrator(Config* cfg);
    void addPotential(InteractionPotential* interactionPot);
    const cx_vec& getMeanField(const int p, const int q);
    virtual void initialize() = 0;
    void computeMeanField(const cx_mat &C);
    virtual cx_vec integrate(const int q, const int r, const cx_mat &C) = 0;
    virtual cx_double integrate(const int p, const int q, const int r, const int s, const cx_mat &C) = 0;
protected:
    vector<InteractionPotential*> potential;

    Config* cfg;

    int nGrid;
    int nOrbitals;

    field<cx_vec> V2;
};

//------------------------------------------------------------------------------
// Inline functions:
//------------------------------------------------------------------------------
inline const cx_vec &MeanFieldIntegrator::getMeanField(const int p, const int q)
{
    return V2(p,q);
}
//------------------------------------------------------------------------------

#endif // MEANFIELDINTEGRATOR_H