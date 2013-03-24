#include "harmonicoscillatoronebody.h"

//------------------------------------------------------------------------------
HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg, const vec &x):
    Potential(cfg, x)
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
        potential(j) = 0.5*w*w*x(j)*x(j);
    }
}
//------------------------------------------------------------------------------
cx_vec HarmonicOscillatorOneBody::evaluate(const cx_vec &psi, double t)
{
    return potential % psi;
}
//------------------------------------------------------------------------------
