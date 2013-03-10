#ifndef COULOMBINTERACTIONNUCLEUS_H
#define COULOMBINTERACTIONNUCLEUS_H

// Local includes
#include <src/Potential/potential.h>

class CoulombInteractionNucleus: public Potential
{
public:
    CoulombInteractionNucleus(Config *cfg);
    virtual cx_vec evaluate(const cx_vec &psi);
};

#endif // COULOMBINTERACTIONNUCLEUS_H
