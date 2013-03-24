#include "potential.h"

//------------------------------------------------------------------------------
Potential::Potential(Config *cfg, const vec &x):
    cfg(cfg),
    x(x)
{
     try{
         nGrid = cfg->lookup("spatialDiscretization.nGrid");
     } catch (const SettingNotFoundException &nfex) {
         cerr << "Potential::Potential(Config *cfg)"
              << "::Error reading from config object." << endl;
     }
}
//------------------------------------------------------------------------------
