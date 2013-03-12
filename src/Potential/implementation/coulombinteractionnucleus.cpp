#include "coulombinteractionnucleus.h"

//------------------------------------------------------------------------------
CoulombInteractionNucleus::CoulombInteractionNucleus(Config *cfg):
    Potential(cfg)
{
    double L;
    double b;
    double Z;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        b = cfg->lookup("oneBodyPotential.coulombInteractionNucleus.b");
        Z = cfg->lookup("oneBodyPotential.coulombInteractionNucleus.Z");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    vec x = linspace(-L, L, nGrid);
    potential = vec(nGrid);

    // Setting the potential
    for(int j=0; j<nGrid; j++){
        potential(j) = - Z/sqrt(x(j)*x(j) + b*b);
    }
}
//------------------------------------------------------------------------------
cx_vec CoulombInteractionNucleus::evaluate(const cx_vec &psi, double t)
{
    return potential % psi ;
}
//------------------------------------------------------------------------------
