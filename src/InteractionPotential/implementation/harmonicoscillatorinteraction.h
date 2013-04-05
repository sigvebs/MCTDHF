#ifndef HARMONICOSCILLATORINTERACTION_H
#define HARMONICOSCILLATORINTERACTION_H

#include <src/InteractionPotential/interactionpotential.h>

class HarmonicOscillatorInteraction: public InteractionPotential
{
public:
    HarmonicOscillatorInteraction(Config* cfg, const Grid &grid);
    virtual double evaluate(uint i, uint j);
protected:
    double epsilon;
};

#endif // HARMONICOSCILLATORINTERACTION_H
