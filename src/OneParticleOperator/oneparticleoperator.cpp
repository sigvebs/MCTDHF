#include "oneparticleoperator.h"
//------------------------------------------------------------------------------
SingleParticleOperator::SingleParticleOperator(Config *cfg, vector<vec> orbitals):
    cfg(cfg),
    orbitals(orbitals)
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
    h = zeros<cx_mat>(nOrbitals, nOrbitals);
    Tspatial = zeros<cx_mat>(nGrid, nSpatialOrbitals);
    Uspatial = zeros<cx_mat>(nGrid, nSpatialOrbitals);
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeNewElements(const cx_mat &C)
{
    this->C = C;
    computeKinetic();
    computePotential();
    computeMatrixElements();
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeMatrixElements()
{
    for(int i=0;i<nOrbitals; i++){
        for(int j=0; j<nOrbitals; j++){
            h(i,j) = 0;
            // Cheking the spin
            if(orbitals[i][0] == orbitals[j][0]){
                h(i,j) += cdot(C.col(i/2), Tspatial.col(j/2));
                h(i,j) += cdot(C.col(i/2), Uspatial.col(j/2));
            }
//            h(j,i) = conj(h(i,j));
        }
    }
#ifdef DEBUG
    cout << "OneParticleOperator::computeMatrixElements()" << endl;
    // Saving kinetic spatial distribution to file.
    //    mat tmp = real(Tspatial);
    //    tmp.save("../DATA/Tspatial.mat", arma_ascii);
//            cout << "h = " << endl << h << endl;
#endif
//            cout << "h = " << endl << h << endl;
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeKinetic()
{
    for(int i=0;i<nSpatialOrbitals; i++){
        for(int j=1; j<nGrid-1; j++){
            Tspatial.col(i)(j) = C(j-1,i) - 2*C(j,i) + C(j+1,i);
        }
        // Endpoints
        Tspatial.col(i)(0) = - 2*C(0,i) + C(1,i);
        Tspatial.col(i)(nGrid-1) = C(nGrid-2,i) - 2*C(nGrid-1,i);
    }
    Tspatial *= -0.5/(dx*dx);
#ifdef DEBUG
    cout << "OneParticleOperator::computeKinetic()" << endl;
    // Saving kinetic spatial distribution to file.
    //    mat tmp = real(Tspatial);
    //    tmp.save("../DATA/Tspatial.mat", arma_ascii);
#endif
    mat tmp = real(Tspatial);
    tmp.save("../DATA/Tspatial.mat", arma_ascii);
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computePotential()
{
    for(int i=0;i<nSpatialOrbitals; i++){
        for(int j=0; j<nGrid; j++){
            Uspatial.col(i)(j) = x[j]*x[j]*C(j,i);
        }
    }
    Uspatial *= 0.5*ww;
#ifdef DEBUG
    cout << "OneParticleOperator::computePotential()" << endl;
    // Saving potential spatial distribution to file.
    //    mat tmp = real(Uspatial);
    //    tmp.save("../DATA/Uspatial.mat", arma_ascii);
#endif
    mat tmp = real(Uspatial);
    tmp.save("../DATA/Uspatial.mat", arma_ascii);
}
//------------------------------------------------------------------------------
const cx_mat &SingleParticleOperator::getH()
{
    return h;
}
//------------------------------------------------------------------------------
const cx_mat &SingleParticleOperator::getHspatial()
{
    return Tspatial + Uspatial;
}
//------------------------------------------------------------------------------
