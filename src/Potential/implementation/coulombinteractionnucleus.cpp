#include "coulombinteractionnucleus.h"

CoulombInteractionNucleus::CoulombInteractionNucleus(Config *cfg):
    Potential(cfg)
{
    double L;
    double b;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        b = cfg->lookup("oneBodyPotential.coulombInteractionNucleus.b");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    vec x = linspace(-L, L, nGrid);
    potential = vec(nGrid);

    // Setting the potential
    for(int j=0; j<nGrid; j++){
        potential(j) = - 2/sqrt(x(j)*x(j) + b*b);
    }
}
//------------------------------------------------------------------------------
cx_vec CoulombInteractionNucleus::evaluate(const cx_vec &psi)
{
    return potential % psi;
}
//------------------------------------------------------------------------------
