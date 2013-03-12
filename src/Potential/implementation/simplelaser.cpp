#include "simplelaser.h"

//------------------------------------------------------------------------------
simpleLaser::simpleLaser(Config *cfg):
    Potential(cfg)
{
    double L;
    double W;
    double e0;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        w = cfg->lookup("oneBodyPotential.simpleLaser.w");
        W = cfg->lookup("oneBodyPotential.harmonicOscillatorBinding.w");
        e0 = cfg->lookup("oneBodyPotential.simpleLaser.e0");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    w *= W;
    potential = linspace(-L, L, nGrid);
    potential *= e0;
}
//------------------------------------------------------------------------------
cx_vec simpleLaser::evaluate(const cx_vec &psi, double t)
{
    return potential % psi * sin(w*t);
}
//------------------------------------------------------------------------------
