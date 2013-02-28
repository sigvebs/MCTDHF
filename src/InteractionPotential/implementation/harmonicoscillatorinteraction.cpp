#include "harmonicoscillatorinteraction.h"


//------------------------------------------------------------------------------
HarmonicOscillatorInteraction::HarmonicOscillatorInteraction(Config *cfg):
    InteractionPotential(cfg)
{
    double L;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        epsilon = cfg->lookup("interactionPotential.interactionPotential.epsilon");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ScreenedCoulombInteraction(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    interactionSpace = zeros(nGrid, nGrid);
    x = linspace(-L, L, nGrid);
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
