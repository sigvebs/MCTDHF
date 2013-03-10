#include "potential.h"

//------------------------------------------------------------------------------
Potential::Potential(Config *cfg):
    cfg(cfg)
{
     try{
         nGrid = cfg->lookup("spatialDiscretization.nGrid");
     } catch (const SettingNotFoundException &nfex) {
         cerr << "Potential::Potential(Config *cfg)"
              << "::Error reading from config object." << endl;
     }
}
//------------------------------------------------------------------------------
