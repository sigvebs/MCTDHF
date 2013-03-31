#include "differentialoperator.h"

//------------------------------------------------------------------------------
DifferentialOperator::DifferentialOperator(Config *cfg, const Grid &grid):
    cfg(cfg),
    grid(grid)
{
     try{
         nGrid = cfg->lookup("spatialDiscretization.nGrid");
     } catch (const SettingNotFoundException &nfex) {
         cerr << "DifferentialOperator::DifferentialOperator(Config *cfg)"
              << "::Error reading from config object." << endl;
     }

    dx = grid.DX;
    dy = grid.DX;
    dz = grid.DZ;
    dxdx = dx*dx;
    dydy = dy*dy;
    dzdz = dz*dz;
}
//------------------------------------------------------------------------------
