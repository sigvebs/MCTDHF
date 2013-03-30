#include "oneparticleoperator.h"
//------------------------------------------------------------------------------
SingleParticleOperator::SingleParticleOperator(Config *cfg, DifferentialOperator *kineticOperator):
    cfg(cfg),
    kineticOperator(kineticOperator)
{
    string filePath;
    try{
        ww = cfg->lookup("oneBodyPotential.harmonicOscillatorBinding.w");
        ww *= ww;
        nOrbitals = cfg->lookup("spatialDiscretization.nSpatialOrbitals");
        nGrid = cfg->lookup("spatialDiscretization.nGrid");
        cfg->lookupValue("systemSettings.filePath", filePath);
    } catch (const SettingNotFoundException &nfex) {
        cerr << "OneParticleOperator::OneParticleOperator(Config *cfg, vector<vec> orbitals)"
             << "::Error reading from config object." << endl;
    }

    h = zeros<cx_mat>(nOrbitals, nOrbitals);
    Tspatial = zeros<cx_mat>(nGrid, nOrbitals);
    Uspatial = zeros<cx_mat>(nGrid, nOrbitals);

    filenameT = filePath + "/T.mat";
    filenameU= filePath + "/U.mat";
}
//------------------------------------------------------------------------------
void SingleParticleOperator::computeNewElements(const cx_mat &C, double t)
{
    computeKinetic(C);
    computePotential(C, t);
    computeMatrixElements(C);
    TU = Tspatial + Uspatial;

#if 1
    for(uint i=0; i<h.n_rows; i++){
        for(uint j=0; j<h.n_rows; j++){
            if(abs(real( h(i,j)) ) < 1e-10)
                h(i,j) = cx_double(0, imag(h(i,j)));

            if(abs(imag( h(i,j)) ) < 1e-10)
                h(i,j) = cx_double(real(h(i,j)), 0);
        }
    }
#endif
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
    Tspatial.save(filenameT);

    // Saving potential spatial distribution to file.
    Uspatial.save(filenameU);
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
    // Cleaning memory
    delete kineticOperator;

    while (!potentials.empty())
     {
       Potential* pot = potentials.back();
       potentials.pop_back();
       delete pot;
     }
}
//------------------------------------------------------------------------------
