#include "oneparticleoperator.h"
//------------------------------------------------------------------------------
SingleParticleOperator::SingleParticleOperator(Config *cfg, DifferentialOperator *kineticOperator):
    cfg(cfg),
    kineticOperator(kineticOperator)
{
    double L;
    try{
        ww = cfg->lookup("oneBodyPotential.harmonicOscillatorBinding.w");
        ww *= ww;
        L = cfg->lookup("spatialDiscretization.latticeRange");
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
    } catch (const SettingNotFoundException &nfex) {
        cerr << "OneParticleOperator::OneParticleOperator(Config *cfg, vector<vec> orbitals)"
             << "::Error reading from config object." << endl;
    }
    x = linspace<cx_vec>(-L, L, nGrid);
    h = zeros<cx_mat>(nOrbitals, nOrbitals);
    Tspatial = zeros<cx_mat>(nGrid, nOrbitals);
    Uspatial = zeros<cx_mat>(nGrid, nOrbitals);
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeNewElements(const cx_mat &C)
{
    // Zeroing out - precaution
    //    Tspatial = zeros<cx_mat>(nGrid, nOrbitals);
    //    Uspatial = zeros<cx_mat>(nGrid, nOrbitals);
    //    h = zeros<cx_mat>(nOrbitals, nOrbitals);
    computeKinetic(C);
    computePotential(C);
    computeMatrixElements(C);
    TU = Tspatial + Uspatial;

    // Saving kinetic spatial distribution to file.
    Tspatial.save("../DATA/Tspatial.mat", arma_ascii);
    // Saving potential spatial distribution to file.
    Uspatial.save("../DATA/Uspatial.mat", arma_ascii);
    //    C.save("../DATA/C.mat", arma_ascii);
}
//------------------------------------------------------------------------------
void SingleParticleOperator::addPotential(Potential *potential)
{
    potentials.push_back(potential);
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeMatrixElements(const cx_mat &C)
{
    h.zeros();

    for(int i=0;i<nOrbitals; i++){
        for(int j=0; j<nOrbitals; j++){
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
    for(int i=0;i<nOrbitals; i++){
        Tspatial.col(i) = kineticOperator->secondDerivative(C.col(i));
    }
    Tspatial *= -0.5;
#ifdef DEBUG
    cout << "OneParticleOperator::computeKinetic()" << endl;
#endif
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computePotential(const cx_mat &C)
{
    Uspatial.zeros();

    for(Potential* potential: potentials){
        for(int i=0;i<nOrbitals; i++){
            Uspatial.col(i) += potential->evaluate(C.col(i));
        }
    }

#ifdef DEBUG
    cout << "OneParticleOperator::computePotential()" << endl;
#endif
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
SingleParticleOperator::~SingleParticleOperator()
{
    delete kineticOperator;
    // TODO: delete all potentials.
}
//------------------------------------------------------------------------------
