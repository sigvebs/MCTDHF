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
}
//------------------------------------------------------------------------------
double HarmonicOscillatorInteraction::evaluate(uint i, uint j)
{
    const vec &r_i = grid.at(i);
    const vec &r_j = grid.at(j);
    return -epsilon*pow(fabs(r_i(0) - r_j(0)), 2);
}
//------------------------------------------------------------------------------
