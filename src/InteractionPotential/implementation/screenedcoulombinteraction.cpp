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
            vec ri = grid.at(i);
            vec rj = grid.at(j);

            interactionSpace(i,j) = lambda/sqrt(pow(rj(0) - ri(0),2) + aa);
        }
    }
    return interactionSpace;
}
//------------------------------------------------------------------------------
