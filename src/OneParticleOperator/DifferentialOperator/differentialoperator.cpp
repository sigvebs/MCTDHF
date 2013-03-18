#include "differentialoperator.h"
//------------------------------------------------------------------------------
DifferentialOperator::DifferentialOperator(Config *cfg):
    cfg(cfg)
{
    double L;
     try{
         L = cfg->lookup("spatialDiscretization.latticeRange");
         nGrid = cfg->lookup("spatialDiscretization.nGrid");
         dx = cfg->lookup("spatialDiscretization.gridSpacing");
     } catch (const SettingNotFoundException &nfex) {
         cerr << "DifferentialOperator::DifferentialOperator(Config *cfg)"
              << "::Error reading from config object." << endl;
     }
    x = linspace<cx_vec>(-L, L - dx, nGrid);
//    dx = real(x(1) - x(0));
    dxdx = dx*dx;
}
//------------------------------------------------------------------------------
