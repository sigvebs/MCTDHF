#include "potential.h"

//------------------------------------------------------------------------------
Potential::Potential(Config *cfg, const Grid &grid):
    cfg(cfg),
    grid(grid)
{
     try{
         nGrid = cfg->lookup("spatialDiscretization.nGrid");
     } catch (const SettingNotFoundException &nfex) {
         cerr << "Potential::Potential(Config *cfg)"
              << "::Error reading from config object." << endl;
         exit(EXIT_FAILURE);
     }
}
//------------------------------------------------------------------------------
