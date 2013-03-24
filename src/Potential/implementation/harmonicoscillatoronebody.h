#ifndef HARMONICOSCILLATORONEBODY_H
#define HARMONICOSCILLATORONEBODY_H

// Local includes
#include <src/Potential/potential.h>

class HarmonicOscillatorOneBody: public Potential
{
public:
    HarmonicOscillatorOneBody(Config *cfg, const vec &x);
    virtual cx_vec evaluate(const cx_vec &psi, double t);
};

#endif // HARMONICOSCILLATORONEBODY_H
