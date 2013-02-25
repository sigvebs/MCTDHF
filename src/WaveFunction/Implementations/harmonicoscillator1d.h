#ifndef HARMONICOSCILLATOR1D_H
#define HARMONICOSCILLATOR1D_H

#include "harmonicoscillator.h"

class HarmonicOscillator1d: public HarmonicOscillator
{
public:
    HarmonicOscillator1d(Config *cfg, vec quantumNumbers);
    virtual mat evaluate(const mat &x);
protected:
    int n;
};

#endif // HARMONICOSCILLATOR1D_H
