#include "harmonicoscillatorinteraction.h"


//------------------------------------------------------------------------------
HarmonicOscillatorInteraction::HarmonicOscillatorInteraction(Config *cfg, const vec &x):
    InteractionPotential(cfg, x)
{
    try{
        epsilon = cfg->lookup("interactionPotential.interactionPotential.epsilon");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ScreenedCoulombInteraction(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    interactionSpace = zeros(nGrid, nGrid);
}
//------------------------------------------------------------------------------
mat HarmonicOscillatorInteraction::computeInteractionSpace()
{
    for(int i=0; i<nGrid; i++){
        for(int j=0; j<nGrid; j++){
            interactionSpace(i,j) = -epsilon*pow(fabs(x(j) - x(i)),2);
        }
    }
    return interactionSpace;
}
//------------------------------------------------------------------------------
