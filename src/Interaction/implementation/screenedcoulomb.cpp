#include "screenedcoulomb.h"
//------------------------------------------------------------------------------
ScreenedCoulomb::ScreenedCoulomb(Config *cfg):
    cfg(cfg)
{
    double L;
    try{
        L = cfg->lookup("spatialDiscretization.latticeRange");
        aa = cfg->lookup("potential.a");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        aa *=aa;
    } catch (const SettingNotFoundException &nfex) {
        cerr << "BasisHarmonicOscillator::BasisHarmonicOscillator(Config *cfg)"
             << "::Error reading from config object." << endl;
    }
    interactionSpace = zeros(nGrid, nGrid);
    x = linspace(-L, L, nGrid);
}
//------------------------------------------------------------------------------
mat ScreenedCoulomb::computeInteractionSpace()
{
    for(int i=0; i<nGrid; i++){
        for(int j=0; j<nGrid; j++){
            interactionSpace(i,j) = 1.0/sqrt(pow(x(j) - x(i),2) + aa);
        }
    }
    return interactionSpace;
}
//------------------------------------------------------------------------------
