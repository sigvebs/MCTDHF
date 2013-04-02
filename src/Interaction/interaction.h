#ifndef INTERACTION_H
#define INTERACTION_H

// Local includes
#include <src/includes/defines.h>
#include <src/includes/lib.h>
#include <src/Interaction/MeanFieldIntegrator/meanfieldintegrator.h>

// Library includes
#include <armadillo>
#include <vector>
#include <unordered_map>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class Interaction
{
public:
    Interaction(Config* cfg, MeanFieldIntegrator* mfIntegrator);
    void computeNewElements(const cx_mat &C);
    const cx_double at(const int p, const int q, const int r, const int s);
    const cx_vec &meanField(const int p, const int q);
    void addPotential(InteractionPotential* interactionPot);
    void updatePositionBasisElements();
    void printInteractionElements();
    ~Interaction();
protected:
    void computeInteractionelements(const cx_mat &C);

    Config* cfg;
    MeanFieldIntegrator* mfIntegrator;

    int nOrbitals;

    // Local storage of data
    unordered_map<int, cx_double> interactionElements;
};

//------------------------------------------------------------------------------
// Inline functions:
//------------------------------------------------------------------------------
inline const cx_vec &Interaction::meanField(const int p, const int q)
{
    return mfIntegrator->getMeanField(p,q);
}
//------------------------------------------------------------------------------
#endif // INTERACTION_H
