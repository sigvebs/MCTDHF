#include "oneparticleoperator.h"
//------------------------------------------------------------------------------
SingleParticleOperator::SingleParticleOperator(Config *cfg, vector<vec> orbitals,
                                               DifferentialOperator *kineticOperator):
    cfg(cfg),
    orbitals(orbitals),
    kineticOperator(kineticOperator)
{
    double L;
    try{
        ww = cfg->lookup("potential.w");
        ww *= ww;
        dx = cfg->lookup("spatialDiscretization.gridSpacing");
        L = cfg->lookup("spatialDiscretization.latticeRange");
        nSpatialOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "OneParticleOperator::OneParticleOperator(Config *cfg, vector<vec> orbitals)"
             << "::Error reading from config object." << endl;
    }
    nOrbitals = orbitals.size();
    x = linspace<cx_vec>(-L, L, nGrid);
    h = zeros<cx_mat>(nSpatialOrbitals, nSpatialOrbitals);
    Tspatial = zeros<cx_mat>(nGrid, nSpatialOrbitals);
    Uspatial = zeros<cx_mat>(nGrid, nSpatialOrbitals);
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeNewElements(const cx_mat &C)
{
    computeKinetic(C);
    computePotential(C);
    computeMatrixElements(C);
    TU = Tspatial + Uspatial;
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeMatrixElements(const cx_mat &C)
{
    for(int i=0;i<nSpatialOrbitals; i++){
        for(int j=0; j<nSpatialOrbitals; j++){
            h(i,j) = 0;
            h(i,j) += cdot(C.col(i), Tspatial.col(j))
                    + cdot(C.col(i), Uspatial.col(j));
        }
    }
#ifdef DEBUG
    cout << "OneParticleOperator::computeMatrixElements()" << endl;
//                cout << "h = " << endl << h << endl;
#endif
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeKinetic(const cx_mat &C)
{
    for(int i=0;i<nSpatialOrbitals; i++){
        Tspatial.col(i) = kineticOperator->secondDerivative(C.col(i));
    }
    Tspatial *= -0.5;
#ifdef DEBUG
    cout << "OneParticleOperator::computeKinetic()" << endl;
    // Saving kinetic spatial distribution to file.
#endif
    Tspatial.save("../DATA/Tspatial.mat", arma_ascii);
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computePotential(const cx_mat &C)
{
    for(int i=0;i<nSpatialOrbitals; i++){
        for(int j=0; j<nGrid; j++){
            Uspatial(j,i) = x(j)*x(j)*C(j,i);
        }
    }
    Uspatial *= 0.5*ww;
#ifdef DEBUG
    cout << "OneParticleOperator::computePotential()" << endl;
    // Saving potential spatial distribution to file.
#endif
    Uspatial.save("../DATA/Uspatial.mat", arma_ascii);
}
//------------------------------------------------------------------------------
const cx_mat &SingleParticleOperator::getH()
{
    return h;
}
//------------------------------------------------------------------------------
const cx_mat &SingleParticleOperator::getHspatial()
{
    return TU;
}
//------------------------------------------------------------------------------
