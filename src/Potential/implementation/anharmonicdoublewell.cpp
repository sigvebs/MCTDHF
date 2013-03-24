#include "anharmonicdoublewell.h"

//------------------------------------------------------------------------------
AnharmonicDoubleWell::AnharmonicDoubleWell(Config *cfg, const vec &x):
    Potential(cfg, x)
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
        potential(j) = pow(x(j) - 0.5*d, 2) * pow(x(j) + 0.5*d, 2);
    }

    potential *= 1.0/(2*d*d);
}
//------------------------------------------------------------------------------
cx_vec AnharmonicDoubleWell::evaluate(const cx_vec &psi, double t)
{
    return potential % psi;
}
//------------------------------------------------------------------------------
