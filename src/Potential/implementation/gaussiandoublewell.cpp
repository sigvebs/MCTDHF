#include "gaussiandoublewell.h"

//------------------------------------------------------------------------------
GaussianDoubleWell::GaussianDoubleWell(Config *cfg, const Grid &grid):
    Potential(cfg, grid)
{
    double V1, V2, Vb, Lx, Ly, Lbx, Lby, a;

    try{
        V1 = cfg->lookup("oneBodyPotential.GaussianDoubleWell.V1");
        V2 = cfg->lookup("oneBodyPotential.GaussianDoubleWell.V2");
        Vb = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Vb");
        Lx = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Lx");
        Ly = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Ly");
        Lbx = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Lbx");
        Lby = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Lby");
        a = cfg->lookup("oneBodyPotential.GaussianDoubleWell.a");

        V2 = V2*V1;
        Lx *= Lx;
        Ly *= Ly;
        Lbx *= Lbx;
        Lby *= Lby;
    } catch (const SettingNotFoundException &nfex) {
        cerr << "GaussianDoubleWell::GaussianDoubleWell(Config *cfg)"
             << "::Error reading from config object." << endl;
        exit(EXIT_FAILURE);
    }
    potential = vec(nGrid);

    // Setting the potential
    for(int j=0; j<nGrid; j++){

        const vec& r = grid.at(j);
        double x = r(0);
        double y = r(1);

        potential(j) =(V1*exp(-pow(x-a,2)/Lx ) + V2*exp(-pow(x+a,2)/Lx) )*exp(-y*y/Ly)
                + Vb*exp(-x*x/Lbx)*exp(-y*y/Lby);
    }
}
//------------------------------------------------------------------------------
cx_vec GaussianDoubleWell::evaluate(const cx_vec &psi, double t)
{
//    cx_vec pot(potential, potential);
//    return pot;
    return potential % psi;
}
//------------------------------------------------------------------------------
