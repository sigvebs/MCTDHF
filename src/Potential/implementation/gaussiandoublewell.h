#ifndef GAUSSIANDOUBLEWELL_H
#define GAUSSIANDOUBLEWELL_H

// Local includes
#include <src/Potential/potential.h>

class GaussianDoubleWell: public Potential
{
public:
    GaussianDoubleWell(Config *cfg,const Grid &grid);
    virtual cx_vec evaluate(const cx_vec &psi, double t);
};

#endif // GAUSSIANDOUBLEWELL_H
