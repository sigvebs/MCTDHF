#include "differentialoperator.h"
//------------------------------------------------------------------------------
DifferentialOperator::DifferentialOperator(Config *cfg):
    cfg(cfg)
{
    double L;
     try{
         dx = cfg->lookup("spatialDiscretization.gridSpacing");
         L = cfg->lookup("spatialDiscretization.latticeRange");
         nGrid = cfg->lookup("spatialDiscretization.nGrid");
     } catch (const SettingNotFoundException &nfex) {
         cerr << "DifferentialOperator::DifferentialOperator(Config *cfg)"
              << "::Error reading from config object." << endl;
     }
    x = linspace<cx_vec>(-L, L, nGrid);
}
//------------------------------------------------------------------------------
