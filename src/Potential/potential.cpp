#include "potential.h"

//------------------------------------------------------------------------------
Potential::Potential(Config *cfg, const Grid &grid):
    grid(grid),
    cfg(cfg)
{
     try{
        dim = cfg->lookup("system.dim");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
     } catch (const SettingNotFoundException &nfex) {
         cerr << "Potential::Potential(Config *cfg)"
              << "::Error reading from config object." << endl;
         exit(EXIT_FAILURE);
     }
}
//------------------------------------------------------------------------------
