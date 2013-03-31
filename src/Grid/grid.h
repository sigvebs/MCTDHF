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
    const vec at(int i) const;
//    double x(int i) const;

    int nGrid = 0;
    int nGridX = 0;
    int nGridY = 0;
    int nGridZ = 0;
    double DX = 0;
    double DY = 0;
    double DZ = 0;
protected:
    mat R;
    int dim;

    void oneDimDiscretization();
    void twoDimDiscretization();
    void threeDimDiscretization();
    Config* cfg;
    string path;
    bool periodicBoundaries;
};

//------------------------------------------------------------------------------
// Inline functions
//------------------------------------------------------------------------------
//inline double Grid::x(int i) const
//{
//    return R(i);
//}
//------------------------------------------------------------------------------
#endif // GRID_H
