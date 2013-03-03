#ifndef ONEPARTICLEOPERATOR_H
#define ONEPARTICLEOPERATOR_H

// Local includes
#include <src/includes/defines.h>
#include <src/includes/lib.h>
#include <src/OneParticleOperator/DifferentialOperator/differentialoperator.h>

// Library includes
#include <armadillo>
#include <vector>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class SingleParticleOperator
{
public:
    SingleParticleOperator(Config* cfg,
                           DifferentialOperator *kineticOperator);
    const cx_mat &getH();
    const cx_mat &getHspatial();
    void computeNewElements(const cx_mat &C);
    ~SingleParticleOperator();
protected:
    void computeMatrixElements(const cx_mat &C);
    void computeKinetic(const cx_mat &C);
    void computePotential(const cx_mat &C);
    Config* cfg;

    // System details
    int nGrid;
    int nOrbitals;
    double dx;
    cx_vec x;

    // Total single particle matrix
    cx_mat h;

    // Kinteic
    DifferentialOperator* kineticOperator;
    cx_mat Tspatial;

    // Potential
    cx_mat Uspatial;
    double ww;

    // Total one-body operator
    cx_mat TU;
};

#endif // ONEPARTICLEOPERATOR_H
