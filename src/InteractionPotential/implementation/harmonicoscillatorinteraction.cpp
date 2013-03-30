#include "harmonicoscillatorinteraction.h"


//------------------------------------------------------------------------------
HarmonicOscillatorInteraction::HarmonicOscillatorInteraction(Config *cfg, const Grid &grid):
    InteractionPotential(cfg, grid)
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
            interactionSpace(i,j) = -epsilon*pow(fabs(grid.x(j) - grid.x(i)),2);
        }
    }
    return interactionSpace;
}
//------------------------------------------------------------------------------
