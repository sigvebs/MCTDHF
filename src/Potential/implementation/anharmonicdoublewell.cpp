#include "anharmonicdoublewell.h"

//------------------------------------------------------------------------------
AnharmonicDoubleWell::AnharmonicDoubleWell(Config *cfg, const Grid &grid):
    Potential(cfg, grid)
{
    double d = 1;
    try{
        d = cfg->lookup("oneBodyPotential.anharmonicDoubleWell.d");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "AnharmonicDoubleWell::AnharmonicDoubleWell(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    potential = vec(nGrid);

    // Setting the potential
    for(int j=0; j<nGrid; j++){
        vec r = grid.at(j);
        potential(j) = pow(r(0) - 0.5*d, 2) * pow(r(0) + 0.5*d, 2);
    }

    potential *= 1.0/(2*d*d);
}
//------------------------------------------------------------------------------
cx_vec AnharmonicDoubleWell::evaluate(const cx_vec &psi, double t)
{
    return potential % psi;
}
//------------------------------------------------------------------------------
