#ifndef GRID_H
#define GRID_H

// Library incldues
#include <armadillo>
#include <vector>
#include <libconfig.h++>

using namespace libconfig;
using namespace std;
using namespace arma;

class Grid
{
public:
    Grid(Config *cfg);
    void createInitalDiscretization();
    void saveGrid();
    void loadGrid();
    double x(int i) const;
    double y(int i) const;
    double z(int i) const;

    vec X;
    vec Y;
    vec Z;
protected:
    int nGrid;
    int dim;
    double DX = 0;
    double DY = 0;
    double DZ = 0;

    Config* cfg;
    string path;
};

//------------------------------------------------------------------------------
// Inline functions
//------------------------------------------------------------------------------
inline double Grid::x(int i) const
{
    return X(i);
}
//------------------------------------------------------------------------------
inline double Grid::y(int i) const
{
    return Y(i);
}
//------------------------------------------------------------------------------
inline double Grid::z(int i) const
{
    return Z(i);
}
//------------------------------------------------------------------------------
#endif // GRID_H
