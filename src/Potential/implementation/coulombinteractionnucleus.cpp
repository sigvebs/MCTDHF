#include "coulombinteractionnucleus.h"

//------------------------------------------------------------------------------
CoulombInteractionNucleus::CoulombInteractionNucleus(Config *cfg, const Grid &grid):
    Potential(cfg, grid)
{
    double b = 1;
    double Z = 2;
    try{
        b = cfg->lookup("oneBodyPotential.coulombInteractionNucleus.b");
        Z = cfg->lookup("oneBodyPotential.coulombInteractionNucleus.Z");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "HarmonicOscillatorOneBody::HarmonicOscillatorOneBody(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    potential = vec(nGrid);

    // Setting the potential
    for(int j=0; j<nGrid; j++){
        potential(j) = - Z/sqrt(grid.x(j)*grid.x(j) + b*b);
    }
}
//------------------------------------------------------------------------------
cx_vec CoulombInteractionNucleus::evaluate(const cx_vec &psi, double t)
{
    return potential % psi ;
}
//------------------------------------------------------------------------------
