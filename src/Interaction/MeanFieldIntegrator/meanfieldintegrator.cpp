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

    for(int p=0; p < nOrbitals; p++){
        for(int q=0; q < nOrbitals; q++){
            V2(p,q) = zeros<cx_vec>(nGrid);
        }
    }
}
//------------------------------------------------------------------------------
void MeanFieldIntegrator::computeMeanField(const cx_mat &C)
{

    for (int q = 0; q < nOrbitals; q++) {

         integrate(q, q, C, V2(q,q));
//        V2(q,q) = integrate(q, q, C);

        for (int r = q+1; r < nOrbitals; r++) {
            integrate(q, r, C, V2(q,r));
//            V2(q,r) = integrate(q, r, C);
            V2(r,q) = conj(V2(q,r));
        }
    }
}
//------------------------------------------------------------------------------
void MeanFieldIntegrator::addPotential(InteractionPotential *interactionPot)
{
    potential.push_back(interactionPot);
}
//------------------------------------------------------------------------------
MeanFieldIntegrator::~MeanFieldIntegrator()
{
    while (!potential.empty())
     {
       InteractionPotential* pot = potential.back();
       potential.pop_back();
       delete pot;
     }
}
//------------------------------------------------------------------------------
