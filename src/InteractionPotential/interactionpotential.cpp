#include "interactionpotential.h"
//------------------------------------------------------------------------------
InteractionPotential::InteractionPotential(Config *cfg, const vec &x):
    cfg(cfg),
    x(x)
{
    try{
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "InteractionPotential(Config *cfg, const vec &x)"
             << "::Error reading from config object." << endl;
    }
}
//------------------------------------------------------------------------------
