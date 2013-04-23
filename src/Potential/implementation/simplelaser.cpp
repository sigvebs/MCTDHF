#include "simplelaser.h"

//------------------------------------------------------------------------------
simpleLaser::simpleLaser(Config *cfg, const Grid &grid):
    Potential(cfg, grid)
{
    double W = 1;
    double e0 = 1;
    int axis = 0;
    try{
        w = cfg->lookup("oneBodyPotential.simpleLaser.w");
        W = cfg->lookup("oneBodyPotential.harmonicOscillatorBinding.w");
        e0 = cfg->lookup("oneBodyPotential.simpleLaser.e0");
        axis = cfg->lookup("oneBodyPotential.simpleLaser.axis");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "simpleLaser::simpleLaser(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    w *= W;
    potential = vec(nGrid);

    // The laser is applied only in the x-direction.
    for(int i=0; i<nGrid; i++){
        vec r = grid.at(i);
        potential(i) = r(axis)*e0;
    }
}
//------------------------------------------------------------------------------
cx_vec simpleLaser::evaluate(const cx_vec &psi, double t)
{
    return potential % psi * sin(w*t);
}
//------------------------------------------------------------------------------
