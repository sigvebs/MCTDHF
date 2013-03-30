#include "simplelaser.h"

//------------------------------------------------------------------------------
simpleLaser::simpleLaser(Config *cfg, const Grid &grid):
    Potential(cfg, grid)
{
    double W = 1;
    double e0 = 1;
    try{
        w = cfg->lookup("oneBodyPotential.simpleLaser.w");
        W = cfg->lookup("oneBodyPotential.harmonicOscillatorBinding.w");
        e0 = cfg->lookup("oneBodyPotential.simpleLaser.e0");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    w *= W;
    potential = grid.X*e0;
}
//------------------------------------------------------------------------------
cx_vec simpleLaser::evaluate(const cx_vec &psi, double t)
{
    return potential % psi * sin(w*t);
}
//------------------------------------------------------------------------------
