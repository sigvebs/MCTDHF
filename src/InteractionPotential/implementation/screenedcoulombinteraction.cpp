#include "screenedcoulombinteraction.h"
//------------------------------------------------------------------------------
ScreenedCoulombInteraction::ScreenedCoulombInteraction(Config *cfg, const Grid &grid):
    InteractionPotential(cfg, grid)
{
    try{
        aa = cfg->lookup("interactionPotential.shieldedCoulombInteraction.a");
        aa *=aa;
        lambda = cfg->lookup("interactionPotential.shieldedCoulombInteraction.lambda");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "ScreenedCoulombInteraction(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
}
//------------------------------------------------------------------------------
double ScreenedCoulombInteraction::evaluate(uint i, uint j)
{
    const vec &r_i = grid.at(i);
    const vec &r_j = grid.at(j);

    double denominator = 0;
    for(int d=0; d<dim; d++){
        denominator += pow(r_i(d) - r_j(d), 2);
    }
    denominator += aa;
    return lambda/sqrt(denominator);
}
//------------------------------------------------------------------------------
