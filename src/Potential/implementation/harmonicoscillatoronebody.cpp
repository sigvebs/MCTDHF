#include "harmonicoscillatoronebody.h"

//------------------------------------------------------------------------------
HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg, const Grid &grid):
    Potential(cfg, grid)
{
    double w = 1;
    try{
        w = cfg->lookup("oneBodyPotential.harmonicOscillatorBinding.w");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    potential = vec(nGrid);

    // Setting the potential
    for(int j=0; j<nGrid; j++){
        const vec& r = grid.at(j);
        double r2 = 0;
        for(int i=0; i<dim; i++){
            r2 += r(i)*r(i);
        }
        potential(j) = 0.5*w*w*r2;
    }
}
//------------------------------------------------------------------------------
cx_vec HarmonicOscillatorOneBody::evaluate(const cx_vec &psi, double t)
{
    return potential % psi;
}
//------------------------------------------------------------------------------
