#ifndef FINITEHARMONICOSCILLATOR_OB_H
#define FINITEHARMONICOSCILLATOR_OB_H

// Local includes
#include <src/Potential/potential.h>

class FiniteHarmonicOscillator_OB: public Potential
{
public:
    FiniteHarmonicOscillator_OB(Config *cfg,const Grid &grid);
    virtual cx_vec evaluate(const cx_vec &psi, double t);
};

#endif // FINITEHARMONICOSCILLATOR_OB_H
