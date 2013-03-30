#include "screenedcoulombinteraction.h"
//------------------------------------------------------------------------------
ScreenedCoulombInteraction::ScreenedCoulombInteraction(Config *cfg, const Grid &grid):
    InteractionPotential(cfg, grid)
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
            interactionSpace(i,j) = lambda/sqrt(pow(grid.x(j) - grid.x(i),2) + aa);
        }
    }
    return interactionSpace;
}
//------------------------------------------------------------------------------
