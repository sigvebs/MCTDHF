#ifndef BASIS_H
#define BASIS_H

// Local includes
#include <src/includes/defines.h>
#include <src/WaveFunction/wavefunction.h>

// Library incldues
#include <armadillo>
#include <vector>
#include <libconfig.h++>


using namespace libconfig;
using namespace std;
using namespace arma;

class Basis
{
public:
    Basis(Config* cfg);
    ~Basis();
    void createBasis();
    void discretizeBasis();
    const vector<vec> &getBasis() const;
    virtual void createInitalDiscretization() = 0;
    const cx_mat &getInitalOrbitals() const;
    const vec &getX() const;
protected:
    void createCartesianBasis();
    void createPolarBasis();

    Config* cfg;

    // Sysmte settings
    int coordinateType;
    int nBasis;
    int dim;

    int nGrid;
    int nSpatialOrbitals;

    Wavefunction *wf;
    vector<vec> states;
    cx_mat C;

    string filnameAxis;

    double dx;
    vec x;
};

#endif // BASIS_H
