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
#ifdef USE_MPI
// disable annoying unused parameter warnings from the MPI library which we don't have any control over
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
// Enable warnings again
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif
using namespace libconfig;
using namespace std;
using namespace arma;

class MeanFieldIntegrator
{
public:
    MeanFieldIntegrator(Config* cfg, const Grid &grid);
    ~MeanFieldIntegrator();
    void addPotential(InteractionPotential* interactionPot);
    const cx_vec& getMeanField(const int p, const int q);
    void computeMeanField(const cx_mat &C);
    virtual void initialize() = 0;
    virtual void integrate(const int q, const int r, const cx_mat &C, cx_vec &V2) = 0;
    virtual cx_double integrate(const int p, const int q, const int r, const int s, const cx_mat &C) = 0;
    void cleanUp();
protected:
    vector<InteractionPotential*> potential;

    Config* cfg;
    const Grid &grid;

    field<cx_vec> V2;
    int nGrid;
    int nOrbitals;

    // MPI
    imat allQR;
    int myRank, nNodes;
    vector<pair<int,int> > myQR;
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
