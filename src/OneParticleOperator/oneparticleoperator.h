#ifndef ONEPARTICLEOPERATOR_H
#define ONEPARTICLEOPERATOR_H

// Local includes
#include <src/includes/defines.h>
#include <src/includes/lib.h>

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
    SingleParticleOperator(Config* cfg, vector<vec> orbitals);
    const cx_mat &getH();
    const cx_mat &getHspatial();
    void computeNewElements(const cx_mat &C);
protected:
    void computeMatrixElements();
    void computeKinetic();
    void computePotential();
    Config* cfg;

    vector<vec> orbitals;
    cx_mat C;

    // System details
    int nGrid;
    int nOrbitals;
    int nSpatialOrbitals;
    double dx;
    cx_vec x;

    // Total single particle matrix
    cx_mat h;

    // Kinteic
    cx_mat Tspatial;

    // Potential
    cx_mat Uspatial;
    double ww;
};

#endif // ONEPARTICLEOPERATOR_H
