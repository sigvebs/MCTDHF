#include "gaussiandoublewell.h"

//------------------------------------------------------------------------------
GaussianDoubleWell::GaussianDoubleWell(Config *cfg, const Grid &grid):
    Potential(cfg, grid)
{
    double V0, VA, VB, VC, Lx, Ly, Lbx, Lby, a;

    try{
        V0 = cfg->lookup("oneBodyPotential.GaussianDoubleWell.V0");
        VA = cfg->lookup("oneBodyPotential.GaussianDoubleWell.VA");
        VB = cfg->lookup("oneBodyPotential.GaussianDoubleWell.VB");
        VC = cfg->lookup("oneBodyPotential.GaussianDoubleWell.VC");
        Lx = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Lx");
        Ly = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Ly");
        Lbx = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Lbx");
        Lby = cfg->lookup("oneBodyPotential.GaussianDoubleWell.Lby");
        a = cfg->lookup("oneBodyPotential.GaussianDoubleWell.a");

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

        potential(j) =(VA*exp(-pow(x-a,2)/Lx ) + VB*exp(-pow(x+a,2)/Lx) )*exp(-y*y/Ly)
                + VC*exp(-x*x/Lbx)*exp(-y*y/Lby);
    }

    potential += V0; // Shift to make the potential positive.
}
//------------------------------------------------------------------------------
cx_vec GaussianDoubleWell::evaluate(const cx_vec &psi, double t)
{
    return potential % psi;
}
//------------------------------------------------------------------------------
