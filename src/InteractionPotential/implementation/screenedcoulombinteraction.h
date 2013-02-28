#ifndef SCREENEDCOULOMBINTERACTION_H
#define SCREENEDCOULOMBINTERACTION_H

#include <src/InteractionPotential/interactionpotential.h>

class ScreenedCoulombInteraction: public InteractionPotential
{
public:
    ScreenedCoulombInteraction(Config* cfg);
    virtual mat computeInteractionSpace();
protected:
    double aa;
};

#endif // SCREENEDCOULOMBINTERACTION_H
