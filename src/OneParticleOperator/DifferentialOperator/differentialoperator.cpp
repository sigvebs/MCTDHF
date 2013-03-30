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
    dx = grid.x(1) - grid.x(0);
    dxdx = dx*dx;
}
//------------------------------------------------------------------------------
