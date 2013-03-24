#include "differentialoperator.h"
//------------------------------------------------------------------------------
DifferentialOperator::DifferentialOperator(Config *cfg, const vec &x):
    cfg(cfg),
    x(x)
{
     try{
         nGrid = cfg->lookup("spatialDiscretization.nGrid");
     } catch (const SettingNotFoundException &nfex) {
         cerr << "DifferentialOperator::DifferentialOperator(Config *cfg)"
              << "::Error reading from config object." << endl;
     }
    dx = x(1) - x(0);
    dxdx = dx*dx;
}
//------------------------------------------------------------------------------
