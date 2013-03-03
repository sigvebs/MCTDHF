#include "meanfieldintegrator.h"

//------------------------------------------------------------------------------
MeanFieldIntegrator::MeanFieldIntegrator(Config *cfg):
    cfg(cfg)
{
    try{
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "MeanFieldIntegrator::MeanFieldIntegrator(Config *cfg)"
             << "::Error reading from config object." << endl;
    }

    V2 = field<cx_vec>(nOrbitals, nOrbitals);
}
//------------------------------------------------------------------------------
void MeanFieldIntegrator::addPotential(InteractionPotential *interactionPot)
{
    potential.push_back(interactionPot);
}
//------------------------------------------------------------------------------
void MeanFieldIntegrator::computeMeanField(const cx_mat &C)
{
    for (int q = 0; q < nOrbitals; q++) {
        for (int r = q; r < nOrbitals; r++) {
            V2(q,r) = integrate(q,r, C);
            V2(r,q) = conj(V2(q,r));
        }
    }
}
//------------------------------------------------------------------------------
