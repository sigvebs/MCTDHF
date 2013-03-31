#ifndef SCREENEDCOULOMBINTERACTION_H
#define SCREENEDCOULOMBINTERACTION_H

#include <src/InteractionPotential/interactionpotential.h>

class ScreenedCoulombInteraction: public InteractionPotential
{
public:
    ScreenedCoulombInteraction(Config* cfg, const Grid &grid);
    virtual mat computeInteractionSpace();
protected:
    double aa;
    double lambda;
};

#endif // SCREENEDCOULOMBINTERACTION_H
