#include "harmonicoscillatoronebody.h"

//------------------------------------------------------------------------------
HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg):
    Potential(cfg)
{
    double L;
    double w;
    double dx;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        w = cfg->lookup("oneBodyPotential.harmonicOscillatorBinding.w");
        dx = cfg->lookup("spatialDiscretization.gridSpacing");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    vec x = linspace(-L, L - dx, nGrid);
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
