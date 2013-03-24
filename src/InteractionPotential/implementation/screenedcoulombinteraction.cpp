#include "screenedcoulombinteraction.h"
//------------------------------------------------------------------------------
ScreenedCoulombInteraction::ScreenedCoulombInteraction(Config *cfg, const vec &x):
    InteractionPotential(cfg, x)
{
    try{
        aa = cfg->lookup("interactionPotential.shieldedCoulombInteraction.a");
        lambda = cfg->lookup("interactionPotential.shieldedCoulombInteraction.lambda");
        aa *=aa;
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ScreenedCoulombInteraction(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    interactionSpace = zeros(nGrid, nGrid);
}
//------------------------------------------------------------------------------
mat ScreenedCoulombInteraction::computeInteractionSpace()
{
    for(int i=0; i<nGrid; i++){
        for(int j=0; j<nGrid; j++){
            interactionSpace(i,j) = lambda/sqrt(pow(x(j) - x(i),2) + aa);
        }
    }
    return interactionSpace;
}
//------------------------------------------------------------------------------
