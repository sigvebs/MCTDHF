#include "finiteharmonicoscillator_ob.h"

//------------------------------------------------------------------------------
FiniteHarmonicOscillator_OB::FiniteHarmonicOscillator_OB(Config *cfg, const Grid &grid):
    Potential(cfg, grid)
{
    double w = 1;
    double r2Cut = 1;
    try{
        w = cfg->lookup("oneBodyPotential.finiteHarmonicOscillator.w");
        r2Cut = cfg->lookup("oneBodyPotential.finiteHarmonicOscillator.rCutoff");
        r2Cut *= r2Cut;

    } catch (const SettingNotFoundException &nfex) {
        cerr << "FiniteHarmonicOscillator_OB::FiniteHarmonicOscillator_OB(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    potential = vec(nGrid);

    // Setting the potential
    for(int j=0; j<nGrid; j++){
        const vec& r = grid.at(j);
        double r2 = 0;
        for(int i=0; i<dim; i++){
            r2 += r(i)*r(i);
        }

        // Checking the cutoff radius
        if(r2 >= r2Cut)
            r2 = r2Cut;

        potential(j) = 0.5*w*w*r2;
    }
}
//------------------------------------------------------------------------------
cx_vec FiniteHarmonicOscillator_OB::evaluate(const cx_vec &psi, double t)
{
    return potential % psi;
}
//------------------------------------------------------------------------------
