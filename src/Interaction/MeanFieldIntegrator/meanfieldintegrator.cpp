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
        V2(q,q) = integrate(q, q, C);

        for (int r = q+1; r < nOrbitals; r++) {
//            cout << "ID = " << omp_get_thread_num() << " " << q << r << endl;
            V2(q,r) = integrate(q, r, C);
            V2(r,q) = conj(V2(q,r));
        }
    }

//    cout << real(V2(0,0)) << endl;
//    exit(1);
}
//------------------------------------------------------------------------------
