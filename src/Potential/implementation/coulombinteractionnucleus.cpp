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
        cerr << "CoulombInteractionNucleus::CoulombInteractionNucleus(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    potential = vec(nGrid);

    // Setting the potential
    for(int j=0; j<nGrid; j++){
        vec r = grid.at(j);

        double r2 = 0;
        for(int i=0; i<dim; i++){
            r2 += r(i)*r(i);
        }
        potential(j) = - Z/sqrt(r2 + b*b);
    }
}
//------------------------------------------------------------------------------
cx_vec CoulombInteractionNucleus::evaluate(const cx_vec &psi, double t)
{
    return potential % psi ;
}
//------------------------------------------------------------------------------
