#include "interactionpotential.h"
//------------------------------------------------------------------------------
InteractionPotential::InteractionPotential(Config *cfg,  const Grid &grid):
    grid(grid),
    cfg(cfg)
{
    try{
        dim = cfg->lookup("system.dim");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "InteractionPotential(Config *cfg,  const Grid &grid)"
             << "::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
}
//------------------------------------------------------------------------------
