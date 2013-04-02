#ifndef SIMPLELASER_H
#define SIMPLELASER_H

// Local includes
#include <src/Potential/potential.h>

class simpleLaser: public Potential
{
public:
    simpleLaser(Config *cfg, const Grid &grid);
    virtual cx_vec evaluate(const cx_vec &psi, double t);
protected:
    double w;
};

#endif // SIMPLELASER_H
