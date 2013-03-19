#include "oneparticleoperator.h"
//------------------------------------------------------------------------------
SingleParticleOperator::SingleParticleOperator(Config *cfg, DifferentialOperator *kineticOperator):
    cfg(cfg),
    kineticOperator(kineticOperator)
{
    double L;
    double dx;
    string filePath;
    try{
        ww = cfg->lookup("oneBodyPotential.harmonicOscillatorBinding.w");
        ww *= ww;
        L = cfg->lookup("spatialDiscretization.latticeRange");
        dx = cfg->lookup("spatialDiscretization.gridSpacing");
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        cfg->lookupValue("systemSettings.filePath", filePath);
    } catch (const SettingNotFoundException &nfex) {
        cerr << "OneParticleOperator::OneParticleOperator(Config *cfg, vector<vec> orbitals)"
             << "::Error reading from config object." << endl;
    }

    x = linspace<cx_vec>(-L, L - dx, nGrid);
    h = zeros<cx_mat>(nOrbitals, nOrbitals);
    Tspatial = zeros<cx_mat>(nGrid, nOrbitals);
    Uspatial = zeros<cx_mat>(nGrid, nOrbitals);

    filenameT = filePath + "T.mat";
    filenameU= filePath + "U.mat";
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeNewElements(const cx_mat &C, double t)
{
    computeKinetic(C);
    computePotential(C, t);
    computeMatrixElements(C);
    TU = Tspatial + Uspatial;
}
//------------------------------------------------------------------------------
void SingleParticleOperator::addPotential(Potential *potential)
{
    potentials.push_back(potential);
}
//------------------------------------------------------------------------------
void SingleParticleOperator::saveOperators()
{
    // Saving kinetic spatial distribution to file.
    Tspatial.save(filenameT, arma_ascii);

    // Saving potential spatial distribution to file.
    Uspatial.save(filenameU, arma_ascii);
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
    Tspatial.zeros();

    for(int i=0;i<nOrbitals; i++){
        Tspatial.col(i) = kineticOperator->secondDerivative(C.col(i));
    }
    Tspatial *= -0.5;
#ifdef DEBUG
    cout << "OneParticleOperator::computeKinetic()" << endl;
#endif
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computePotential(const cx_mat &C, double t)
{
    Uspatial.zeros();

    for(Potential* potential: potentials){
        for(int i=0;i<nOrbitals; i++){
            Uspatial.col(i) += potential->evaluate(C.col(i), t);
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
