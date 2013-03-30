#ifndef BASIS_H
#define BASIS_H

// Local includes
#include <src/includes/defines.h>
#include <src/Grid/grid.h>
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
    const vector<vec> &getBasis() const;
    const field<cx_mat> &getOrbitals() const;
    void setGrid(Grid *grid);
    void loadOrbitals();
    virtual void createInitalDiscretization() = 0;
protected:
    void createCartesianBasis();
    void createPolarBasis();
    void SVD(cx_mat &D);

    vector<vec> states;
    field<cx_mat> C;

    // Sysmte settings
    int coordinateType;
    int nBasis;
    int dim;
    int nGrid;
    int nSpatialOrbitals;
    Grid *grid;

    string filnameAxis;
    Config* cfg;
};

#endif // BASIS_H
